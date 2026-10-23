//! The inverted index as flat sorted arrays, read in place from memory-mapped files.
//!
//! `ProteomeIndex::finalize` writes four files into `<index>/postings/`, all little-endian:
//!
//! - `shards.u64`: `n_shards + 1` positions in `kmers.u64`. Shard `s` holds the k-mers
//!   `shards[s]..shards[s + 1]`. A k-mer's shard is its hash's low bits.
//! - `kmers.u64`: every distinct k-mer hash, sorted within its shard.
//! - `offsets.u64`: `n_kmers + 1` positions in `targets.u32`. K-mer `i` is found in the
//!   targets `offsets[i]..offsets[i + 1]`, its posting list.
//! - `targets.u32`: target indexes into the target list, each posting list in target order.
//!
//! A search maps the files rather than copying them into hash maps. Opening validates
//! the shard and offset tables. K-mer and target pages are read when lookups touch them;
//! searches running at the same time on one node share mapped pages.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;

use crate::errors::{IndexError, IndexResult};

#[cfg(not(target_endian = "little"))]
compile_error!("the posting files are little-endian and are read in place");

/// Directory inside the index that holds the posting files.
const DIR: &str = "postings";
/// Where `PostingWriter` builds the files before swapping them in.
const TMP_DIR: &str = "postings.tmp";
/// Where the previous files go while the new ones are swapped in.
const OLD_DIR: &str = "postings.old";

const SHARDS_FILE: &str = "shards.u64";
const KMERS_FILE: &str = "kmers.u64";
const OFFSETS_FILE: &str = "offsets.u64";
const TARGETS_FILE: &str = "targets.u32";

/// One posting file mapped into memory. An empty file is not mapped, because mapping
/// zero bytes fails on some platforms.
#[derive(Debug)]
struct Section(Option<Mmap>);

impl Section {
    fn open(dir: &Path, name: &str) -> IndexResult<Self> {
        let path = dir.join(name);
        let file = File::open(&path).map_err(|e| {
            IndexError::CorruptIndex(format!(
                "posting file {} cannot be opened ({e}); rebuild the index with `kmerseek index`",
                path.display()
            ))
        })?;
        if file.metadata()?.len() == 0 {
            return Ok(Self(None));
        }
        // SAFETY: posting files are never written in place. `PostingWriter` builds new ones
        // in another directory and renames it over this one, which leaves these bytes
        // unchanged for as long as they stay mapped.
        let map = unsafe { Mmap::map(&file)? };
        Ok(Self(Some(map)))
    }

    /// The file as a slice of `T`. Never fails after `PostingIndex::open` checked it:
    /// mappings start on a page boundary, so only a length that is not a multiple of
    /// `size_of::<T>()` can be refused.
    fn slice<T: bytemuck::Pod>(&self) -> &[T] {
        // An empty `&[u8]` has no alignment to speak of, so it cannot be cast; return an
        // empty `&[T]` directly.
        self.0.as_deref().map_or(&[], bytemuck::cast_slice)
    }

    fn check<T: bytemuck::Pod>(&self, name: &str) -> IndexResult<()> {
        let Some(bytes) = self.0.as_deref() else { return Ok(()) };
        bytemuck::try_cast_slice::<u8, T>(bytes).map(|_| ()).map_err(|e| {
            IndexError::CorruptIndex(format!("posting file {name} is not a whole array ({e})"))
        })
    }
}

/// The inverted index of a finalized `ProteomeIndex`: for each k-mer hash, the targets
/// that contain it.
#[derive(Debug)]
pub struct PostingIndex {
    shards: Section,
    kmers: Section,
    offsets: Section,
    targets: Section,
}

impl PostingIndex {
    /// Whether the index at `index_dir` has posting files.
    pub fn exists(index_dir: &Path) -> bool {
        index_dir.join(DIR).join(SHARDS_FILE).exists()
    }

    /// Map posting files and validate shard and posting-list boundaries.
    /// This reads the full shard and offset tables, but not all k-mers or targets.
    /// Invalid interior offsets must fail before a search can silently miss hits.
    pub fn open(index_dir: &Path) -> IndexResult<Self> {
        let dir = index_dir.join(DIR);
        let index = Self {
            shards: Section::open(&dir, SHARDS_FILE)?,
            kmers: Section::open(&dir, KMERS_FILE)?,
            offsets: Section::open(&dir, OFFSETS_FILE)?,
            targets: Section::open(&dir, TARGETS_FILE)?,
        };
        index.check()?;
        Ok(index)
    }

    /// An index with no k-mers, for a searcher built without a database.
    pub fn empty() -> Self {
        Self {
            shards: Section(None),
            kmers: Section(None),
            offsets: Section(None),
            targets: Section(None),
        }
    }

    fn check(&self) -> IndexResult<()> {
        self.shards.check::<u64>(SHARDS_FILE)?;
        self.kmers.check::<u64>(KMERS_FILE)?;
        self.offsets.check::<u64>(OFFSETS_FILE)?;
        self.targets.check::<u32>(TARGETS_FILE)?;
        let (shards, offsets) = (self.shards.slice::<u64>(), self.offsets.slice::<u64>());
        let n_shards = shards.len().saturating_sub(1);
        let corrupt =
            |what: String| Err(IndexError::CorruptIndex(format!("posting files: {what}")));
        if !n_shards.is_power_of_two() {
            return corrupt(format!("{n_shards} shards, not a power of two"));
        }
        if shards[0] != 0 || shards[n_shards] != self.len() as u64 {
            return corrupt(format!(
                "shards end at {}, kmers.u64 holds {}",
                shards[n_shards],
                self.len()
            ));
        }
        if offsets.len() != self.len() + 1 {
            return corrupt(format!("{} offsets for {} k-mers", offsets.len(), self.len()));
        }
        if offsets[0] != 0 || offsets[self.len()] != self.n_postings() as u64 {
            return corrupt(format!(
                "offsets end at {}, targets.u32 holds {}",
                offsets[self.len()],
                self.n_postings()
            ));
        }
        // Checking endpoints alone permits corrupt interior boundaries to panic while
        // slicing a shard or to turn a real posting list into a missing lookup.
        if shards.windows(2).any(|w| w[0] > w[1] || w[1] > self.len() as u64) {
            return corrupt("shard boundaries are not monotone and within the k-mer array".into());
        }
        if offsets.windows(2).any(|w| w[0] > w[1] || w[1] > self.n_postings() as u64) {
            return corrupt("posting offsets are not monotone and within the target array".into());
        }
        Ok(())
    }

    /// Distinct k-mer hashes in the index.
    pub fn len(&self) -> usize {
        self.kmers.slice::<u64>().len()
    }

    /// Whether the index holds no k-mers.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Total postings: the sum, over every k-mer, of how many targets contain it.
    pub fn n_postings(&self) -> usize {
        self.targets.slice::<u32>().len()
    }

    /// Number of shards, 0 for `empty()`.
    pub fn n_shards(&self) -> usize {
        self.shards.slice::<u64>().len().saturating_sub(1)
    }

    /// The targets that contain `hash`, in target order; `None` if no target does.
    pub fn get(&self, hash: u64) -> Option<&[u32]> {
        let n_shards = self.n_shards();
        if n_shards == 0 {
            return None;
        }
        let shard = (hash as usize) & (n_shards - 1);
        let (start, end) = self.shard_range(shard);
        let i = start + self.kmers.slice::<u64>()[start..end].binary_search(&hash).ok()?;
        self.posting_list(i)
    }

    /// How many targets contain `hash`: the length of its posting list.
    pub fn frequency(&self, hash: u64) -> Option<usize> {
        self.get(hash).map(<[u32]>::len)
    }

    /// Every `(hash, posting list)` of one shard, in hash order.
    pub fn shard(&self, shard: usize) -> impl Iterator<Item = (u64, &[u32])> + '_ {
        let (start, end) = if shard < self.n_shards() { self.shard_range(shard) } else { (0, 0) };
        let kmers = self.kmers.slice::<u64>();
        (start..end).filter_map(move |i| Some((kmers[i], self.posting_list(i)?)))
    }

    fn shard_range(&self, shard: usize) -> (usize, usize) {
        let shards = self.shards.slice::<u64>();
        (shards[shard] as usize, shards[shard + 1] as usize)
    }

    fn posting_list(&self, i: usize) -> Option<&[u32]> {
        let offsets = self.offsets.slice::<u64>();
        self.targets.slice::<u32>().get(offsets[i] as usize..offsets[i + 1] as usize)
    }
}

/// Writes the posting files one shard at a time, so memory holds one shard, never the
/// whole index.
pub struct PostingWriter {
    index_dir: PathBuf,
    shard_starts: Vec<u64>,
    kmers: BufWriter<File>,
    offsets: BufWriter<File>,
    targets: BufWriter<File>,
    n_kmers: u64,
    n_postings: u64,
}

impl PostingWriter {
    /// Start new posting files for the index at `index_dir`. The current files, if any,
    /// stay readable until `finish`.
    pub fn create(index_dir: &Path) -> IndexResult<Self> {
        let tmp = index_dir.join(TMP_DIR);
        if tmp.exists() {
            fs::remove_dir_all(&tmp)?;
        }
        fs::create_dir(&tmp)?;
        let open = |name: &str| -> IndexResult<BufWriter<File>> {
            Ok(BufWriter::with_capacity(1 << 20, File::create(tmp.join(name))?))
        };
        let mut writer = Self {
            index_dir: index_dir.to_path_buf(),
            shard_starts: vec![0],
            kmers: open(KMERS_FILE)?,
            offsets: open(OFFSETS_FILE)?,
            targets: open(TARGETS_FILE)?,
            n_kmers: 0,
            n_postings: 0,
        };
        writer.offsets.write_all(&0u64.to_le_bytes())?;
        Ok(writer)
    }

    /// Append the next shard: `(hash, posting list)` in hash order.
    pub fn push_shard(&mut self, groups: &[(u64, Vec<u32>)]) -> IndexResult<()> {
        for (hash, targets) in groups {
            self.kmers.write_all(&hash.to_le_bytes())?;
            for target in targets {
                self.targets.write_all(&target.to_le_bytes())?;
            }
            self.n_postings += targets.len() as u64;
            self.offsets.write_all(&self.n_postings.to_le_bytes())?;
        }
        self.n_kmers += groups.len() as u64;
        self.shard_starts.push(self.n_kmers);
        Ok(())
    }

    /// Write the shard table, flush everything to disk, and swap the new files in.
    pub fn finish(self) -> IndexResult<()> {
        let tmp = self.index_dir.join(TMP_DIR);
        let mut shards = BufWriter::new(File::create(tmp.join(SHARDS_FILE))?);
        for start in &self.shard_starts {
            shards.write_all(&start.to_le_bytes())?;
        }
        for writer in [shards, self.kmers, self.offsets, self.targets] {
            writer.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        }

        let (dir, old) = (self.index_dir.join(DIR), self.index_dir.join(OLD_DIR));
        if old.exists() {
            fs::remove_dir_all(&old)?;
        }
        if dir.exists() {
            fs::rename(&dir, &old)?;
        }
        fs::rename(&tmp, &dir)?;
        if old.exists() {
            fs::remove_dir_all(&old)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    /// Two shards: hashes 2 and 4 are even (shard 0), 3 and 7 odd (shard 1).
    fn write_two_shards(dir: &Path) -> IndexResult<()> {
        let mut writer = PostingWriter::create(dir)?;
        writer.push_shard(&[(2, vec![0, 5]), (4, vec![1])])?;
        writer.push_shard(&[(3, vec![2, 3, 4]), (7, vec![0])])?;
        writer.finish()
    }

    #[test]
    fn lookups_find_each_posting_list_and_nothing_else() -> IndexResult<()> {
        let dir = tempdir()?;
        write_two_shards(dir.path())?;
        let index = PostingIndex::open(dir.path())?;
        assert_eq!((index.len(), index.n_postings(), index.n_shards()), (4, 7, 2));
        assert_eq!(index.get(2), Some(&[0u32, 5][..]));
        assert_eq!(index.get(4), Some(&[1u32][..]));
        assert_eq!(index.get(3), Some(&[2u32, 3, 4][..]));
        assert_eq!(index.get(7), Some(&[0u32][..]));
        assert_eq!(index.get(6), None, "even, so shard 0, which has no 6");
        assert_eq!(index.get(5), None, "odd, so shard 1, which has no 5");
        assert_eq!(index.frequency(3), Some(3));
        let shard_1: Vec<(u64, Vec<u32>)> = index.shard(1).map(|(h, t)| (h, t.to_vec())).collect();
        assert_eq!(shard_1, vec![(3, vec![2, 3, 4]), (7, vec![0])]);
        Ok(())
    }

    /// Rewriting the files while the old ones are mapped leaves the old mapping intact,
    /// which is what lets `finalize` read the old shards while writing the new ones.
    #[test]
    fn rewriting_swaps_in_new_files_and_keeps_old_mapping_readable() -> IndexResult<()> {
        let dir = tempdir()?;
        write_two_shards(dir.path())?;
        let old = PostingIndex::open(dir.path())?;
        let mut writer = PostingWriter::create(dir.path())?;
        writer.push_shard(&[(8, vec![9])])?;
        writer.push_shard(&[])?;
        writer.finish()?;
        assert_eq!(old.get(3), Some(&[2u32, 3, 4][..]));
        let new = PostingIndex::open(dir.path())?;
        assert_eq!((new.len(), new.get(8), new.get(3)), (1, Some(&[9u32][..]), None));
        assert!(!dir.path().join(TMP_DIR).exists() && !dir.path().join(OLD_DIR).exists());
        Ok(())
    }

    #[test]
    fn an_index_with_no_kmers_opens_and_finds_nothing() -> IndexResult<()> {
        let dir = tempdir()?;
        let mut writer = PostingWriter::create(dir.path())?;
        writer.push_shard(&[])?;
        writer.finish()?;
        let index = PostingIndex::open(dir.path())?;
        assert_eq!((index.len(), index.n_postings(), index.get(0)), (0, 0, None));
        assert_eq!(PostingIndex::empty().get(0), None);
        Ok(())
    }

    #[test]
    fn a_missing_or_cut_file_is_reported_as_corrupt() -> IndexResult<()> {
        let dir = tempdir()?;
        write_two_shards(dir.path())?;
        let targets = dir.path().join(DIR).join(TARGETS_FILE);
        fs::write(&targets, [0u8; 4])?;
        let message = PostingIndex::open(dir.path()).unwrap_err().to_string();
        assert_eq!(message, "Corrupt index: posting files: offsets end at 7, targets.u32 holds 1");
        fs::write(&targets, [0u8; 5])?;
        let message = PostingIndex::open(dir.path()).unwrap_err().to_string();
        assert!(
            message.starts_with("Corrupt index: posting file targets.u32 is not a whole array"),
            "{message}"
        );
        fs::remove_file(&targets)?;
        let message = PostingIndex::open(dir.path()).unwrap_err().to_string();
        assert!(message.contains("targets.u32 cannot be opened"), "{message}");
        Ok(())
    }
    #[test]
    fn invalid_interior_shard_is_reported_before_lookup() -> IndexResult<()> {
        let dir = tempdir()?;
        write_two_shards(dir.path())?;
        let bytes: Vec<_> = [0u64, 99, 4].into_iter().flat_map(u64::to_le_bytes).collect();
        fs::write(dir.path().join(DIR).join(SHARDS_FILE), bytes)?;
        let error = PostingIndex::open(dir.path()).unwrap_err().to_string();
        assert!(error.contains("shard boundaries"), "{error}");
        Ok(())
    }

    #[test]
    fn invalid_interior_offsets_are_not_silently_missing_hits() -> IndexResult<()> {
        let dir = tempdir()?;
        write_two_shards(dir.path())?;
        for offsets in [[0u64, 99, 3, 6, 7], [0, 2, 1, 6, 7]] {
            let bytes: Vec<_> = offsets.into_iter().flat_map(u64::to_le_bytes).collect();
            fs::write(dir.path().join(DIR).join(OFFSETS_FILE), bytes)?;
            let error = PostingIndex::open(dir.path()).unwrap_err().to_string();
            assert!(error.contains("posting offsets"), "{error}");
        }
        Ok(())
    }
}
