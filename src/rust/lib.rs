use pyo3::prelude::*;
use pyo3::types::PyType;
use sourmash::encodings::{self, HashFunctions};
use std::path::PathBuf;

pub mod aminoacid;
pub mod index;
pub mod uniprot;

#[cfg(test)]
mod tests;

// Re-export main types for easier access
pub use aminoacid::AminoAcidAmbiguity;
pub use index::ProteomeIndex;
pub use uniprot::{ProteinFeature, UniProtEntry};

#[pyclass]
#[derive(Clone, Copy)]
pub enum PyProteinEncoding {
    Raw,
    Dayhoff,
    HP,
}

#[pymethods]
impl PyProteinEncoding {
    #[classmethod]
    fn raw(_cls: Py<PyType>) -> Self {
        Self::Raw
    }

    #[classmethod]
    fn dayhoff(_cls: Py<PyType>) -> Self {
        Self::Dayhoff
    }

    #[classmethod]
    fn hp(_cls: Py<PyType>) -> Self {
        Self::HP
    }
}

impl From<PyProteinEncoding> for String {
    fn from(py_encoding: PyProteinEncoding) -> Self {
        match py_encoding {
            PyProteinEncoding::Raw => "protein".to_string(),
            PyProteinEncoding::Dayhoff => "dayhoff".to_string(),
            PyProteinEncoding::HP => "hp".to_string(),
        }
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyclass]
pub struct PyProteomeIndex {
    index: crate::index::ProteomeIndex,
}

#[pymethods]
impl PyProteomeIndex {
    #[new]
    fn new(
        ksize: u32,
        scaled: u32,
        moltype: &str,
        encoding_type: &str,
        db_path: &str,
    ) -> PyResult<Self> {
        // Convert scaled to u64
        let scaled = scaled as u64;

        // Convert moltype string to HashFunctions enum
        let hash_function = match moltype {
            "protein" => HashFunctions::Murmur64Protein,
            "dayhoff" => HashFunctions::Murmur64Dayhoff,
            "hp" => HashFunctions::Murmur64Hp,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Invalid moltype",
                ))
            }
        };

        // Convert encoding_type to encoding function
        let encoding_fn: fn(u8) -> u8 = match encoding_type {
            "raw" => |x| x,
            "dayhoff" => sourmash::encodings::aa_to_dayhoff,
            "hp" => sourmash::encodings::aa_to_hp,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Invalid encoding type",
                ))
            }
        };

        let index = crate::index::ProteomeIndex::new(
            PathBuf::from(db_path),
            ksize,
            scaled,
            hash_function,
            encoding_fn,
        )
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        Ok(PyProteomeIndex { index })
    }

    fn process_uniprot_xml(&self, path: &str) -> PyResult<()> {
        self.index
            .process_uniprot_xml(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))
    }

    fn process_protein_files(&self, files: Vec<String>) -> PyResult<()> {
        let paths: Vec<PathBuf> = files.into_iter().map(PathBuf::from).collect();
        self.index
            .process_protein_files(&paths)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))
    }

    fn get_mins(&self) -> Vec<u64> {
        self.index.get_hashes()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn kmerseek(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<PyProteinEncoding>()?;
    m.add_class::<PyProteomeIndex>()?;
    Ok(())
}
