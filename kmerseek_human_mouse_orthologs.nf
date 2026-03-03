#!/usr/bin/env nextflow

/*
 * Nextflow pipeline to evaluate kmerseek ortholog detection using human-mouse orthologs
 *
 * This pipeline:
 * 1. Downloads the JAX ortholog mapping file
 * 2. Indexes human and mouse GENCODE protein sequences
 * 3. Searches human proteins against mouse proteins for k=15-30
 * 4. Evaluates ortholog detection accuracy at each k-size
 *
 * Usage:
 *   nextflow run kmerseek_human_mouse_orthologs.nf
 */

// FASTA files for human and mouse protein sequences
params.human_fasta = "${System.getProperty('user.home')}/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa"
params.mouse_fasta = "${System.getProperty('user.home')}/data/gencode/mouse/m38/gencode.vM38.pc_translations.canonical.fa"
params.outdir = "${System.getProperty('user.home')}/data/gencode/results-human-mouse-orthologs"
params.kmerseek = "${System.getProperty('user.home')}/code/kmerseek/target/release/kmerseek-rust"

// Ortholog mapping URL
params.ortholog_url = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"


// Minimum containment threshold — 0.0 keeps all hits (large CSVs; polars handles them)
params.threshold = 0.0

process buildRelease {
    output:
    path 'kmerseek-rust'

    script:
    def kmerseek_dir = "${System.getProperty('user.home')}/code/kmerseek"
    """
    WORK_DIR=\$PWD
    cd ${kmerseek_dir}
    cargo build --release
    cp ${kmerseek_dir}/target/release/kmerseek-rust \$WORK_DIR/kmerseek-rust
    """
}

process downloadOrthologMapping {
    publishDir params.outdir, mode: 'copy'

    output:
    path 'HOM_MouseHumanSequence.rpt'

    script:
    """
    curl -o HOM_MouseHumanSequence.rpt "${params.ortholog_url}"
    """
}

process parseOrthologMapping {
    publishDir params.outdir, mode: 'copy'

    input:
    path ortholog_file

    output:
    path 'ortholog_pairs.tsv'
    path 'ortholog_stats.txt'

    script:
    """
    #!/usr/bin/env python3
    import csv
    from collections import defaultdict

    # Parse ortholog file and extract human-mouse gene symbol pairs
    # Group by DB Class Key to identify ortholog groups
    ortholog_groups = defaultdict(lambda: {'human': set(), 'mouse': set()})

    with open('${ortholog_file}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            db_class_key = row['DB Class Key']
            organism = row['Common Organism Name']
            symbol = row['Symbol']

            if organism == 'human':
                # Human gene symbols are uppercase in GENCODE
                ortholog_groups[db_class_key]['human'].add(symbol.upper())
            elif organism == 'mouse, laboratory':
                # Mouse gene symbols are capitalized (first letter uppercase) in GENCODE
                # But the JAX file already has proper capitalization
                ortholog_groups[db_class_key]['mouse'].add(symbol)

    # Create pairs: for each human gene, list all mouse orthologs
    human_to_mouse = defaultdict(set)
    mouse_to_human = defaultdict(set)

    for group in ortholog_groups.values():
        for human_gene in group['human']:
            for mouse_gene in group['mouse']:
                human_to_mouse[human_gene].add(mouse_gene)
                mouse_to_human[mouse_gene].add(human_gene)

    # Write pairs file (human_gene, mouse_gene)
    with open('ortholog_pairs.tsv', 'w') as f:
        f.write('human_gene\\tmouse_gene\\n')
        for human_gene, mouse_genes in sorted(human_to_mouse.items()):
            for mouse_gene in sorted(mouse_genes):
                f.write(f'{human_gene}\\t{mouse_gene}\\n')

    # Write stats
    with open('ortholog_stats.txt', 'w') as f:
        f.write(f'Number of ortholog groups: {len(ortholog_groups)}\\n')
        f.write(f'Number of human genes with mouse orthologs: {len(human_to_mouse)}\\n')
        f.write(f'Number of mouse genes with human orthologs: {len(mouse_to_human)}\\n')
        f.write(f'Total human-mouse pairs: {sum(len(v) for v in human_to_mouse.values())}\\n')

        # Check for one-to-many and many-to-many relationships
        one_to_one = sum(1 for v in human_to_mouse.values() if len(v) == 1)
        one_to_many = sum(1 for v in human_to_mouse.values() if len(v) > 1)
        f.write(f'Human genes with exactly one mouse ortholog: {one_to_one}\\n')
        f.write(f'Human genes with multiple mouse orthologs: {one_to_many}\\n')
    """
}

process decompressFasta {
    tag "${species}_${fasta.name}"

    input:
    tuple val(species), path(fasta)

    output:
    tuple val(species), path("${fasta.baseName}")

    script:
    """
    gunzip -c ${fasta} > ${fasta.baseName}
    """
}

process indexDatabase {
    tag "${species}_hp_k${ksize}"
    publishDir "${params.outdir}/indices", mode: 'copy', pattern: '*.rocksdb', type: 'dir'
    publishDir params.outdir, mode: 'copy', pattern: '*.index.log'

    input:
    path kmerseek
    tuple val(species), path(fasta), val(ksize)

    output:
    tuple val(species), val(ksize), path("${fasta}.hp.k${ksize}.scaled1.kmerseek.rocksdb", type: 'dir')
    path "${fasta}.hp.k${ksize}.scaled1.kmerseek.index.log"

    script:
    def log_file = "${fasta}.hp.k${ksize}.scaled1.kmerseek.index.log"
    """
    echo "=== Indexing: ${species} hp k=${ksize} ===" | tee ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} index \\
        --encoding hp \\
        --ksize ${ksize} \\
        --scaled 1 \\
        --input ${fasta} \\
        2>&1 | tee -a ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    """
}

process searchHumanVsMouse {
    tag "hp_k${ksize}"
    publishDir params.outdir, mode: 'copy', pattern: '*.csv'
    publishDir params.outdir, mode: 'copy', pattern: '*.search.log'

    input:
    path kmerseek
    tuple val(ksize), path(human_fasta), path(mouse_index)

    output:
    tuple val(ksize), path("human_vs_mouse.hp.k${ksize}.results.csv")
    path "human_vs_mouse.hp.k${ksize}.search.log"

    script:
    def output_csv = "human_vs_mouse.hp.k${ksize}.results.csv"
    def log_file = "human_vs_mouse.hp.k${ksize}.search.log"
    """
    echo "=== Searching: human vs mouse hp k=${ksize} ===" | tee ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} search \\
        --encoding hp \\
        --ksize ${ksize} \\
        --query ${human_fasta} \\
        --target ${mouse_index} \\
        > ${output_csv} 2>> ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "Results: \$(wc -l < ${output_csv}) rows" | tee -a ${log_file}
    """
}

process evaluateOrthologs {
    tag "hp_k${ksize}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(ksize), path(results_csv), path(ortholog_pairs)

    output:
    tuple val(ksize), path("ortholog_evaluation.hp.k${ksize}.tsv")
    path "ortholog_evaluation.hp.k${ksize}.summary.txt"
    path "ortholog_evaluation.hp.k${ksize}.roc_data.tsv"
    path "ortholog_evaluation.hp.k${ksize}.mht.csv"
    path "metrics_*.hp.k${ksize}.png"

    script:
    """
    #!/usr/bin/env python3
    import polars as pl
    import numpy as np
    import json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ksize = ${ksize}

    # ── Load ortholog ground truth ─────────────────────────────────────────────
    ortho = (
        pl.read_csv('${ortholog_pairs}', separator='\\t')
        .with_columns(pl.col('human_gene').str.to_uppercase(), pl.lit(True).alias('is_ortholog'))
    )

    # ── Lazy scan results + gene-symbol extraction + ortholog join ─────────────
    df = (
        pl.scan_csv('${results_csv}')
        .with_columns([
            pl.col('query_name').str.split('|').list.get(6).alias('human_gene'),
            pl.col('target_name').str.split('|').list.get(6).alias('mouse_gene'),
        ])
        .with_columns(pl.col('human_gene').str.to_uppercase())
        .join(
            ortho.lazy().select(['human_gene', 'mouse_gene', 'is_ortholog']),
            on=['human_gene', 'mouse_gene'], how='left',
        )
        .with_columns(pl.col('is_ortholog').fill_null(False))
        .collect()
    )

    n_total = len(df)
    n_orth  = int(df['is_ortholog'].sum())
    n_non   = n_total - n_orth
    print(f"Loaded {n_total:,} rows  |  {n_orth:,} orthologs  |  {n_non:,} non-orthologs")

    # ── Metric columns ─────────────────────────────────────────────────────────
    METRIC_COLS = [
        'n_intersecting_hashes', 'jaccard', 'containment', 'query_tfidf',
        'mean_matched_kmer_freq', 'sum_matched_kmer_freq',
        'expected_shared_kmers', 'enrichment', 'prob_overlap',
    ]
    metrics = [c for c in METRIC_COLS if c in df.columns]

    # ── Histograms (ortholog vs non-ortholog) ──────────────────────────────────
    orth_df = df.filter(pl.col('is_ortholog'))
    nort_df = df.filter(~pl.col('is_ortholog'))

    for metric in metrics:
        ov  = orth_df[metric].drop_nulls().to_numpy()
        nov = nort_df[metric].drop_nulls().to_numpy()
        if len(ov) == 0 and len(nov) == 0:
            continue
        all_vals = np.concatenate([ov, nov])
        lo, hi   = np.nanpercentile(all_vals, [0.5, 99.5])
        bins     = np.linspace(lo, hi, 60)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(nov, bins=bins, alpha=0.6, density=True,
                label=f'Non-ortholog (n={len(nov):,})', color='steelblue')
        ax.hist(ov,  bins=bins, alpha=0.6, density=True,
                label=f'Ortholog (n={len(ov):,})', color='tomato')
        ax.set_xlabel(metric)
        ax.set_ylabel('Density')
        ax.set_title(f'{metric}  —  HP k={ksize}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(f'metrics_{metric}.hp.k{ksize}.png', dpi=150)
        plt.close()

    # ── MHT corrections on prob_overlap ───────────────────────────────────────
    pvals = df['prob_overlap'].fill_null(1.0).to_numpy()
    n = len(pvals)

    def bh_adj(p):
        idx = np.argsort(p); sp = p[idx]; m = len(p)
        adj = np.minimum(1, sp * m / np.arange(1, m + 1))
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        out = np.empty(m); out[idx] = adj; return out

    def by_adj(p):
        idx = np.argsort(p); sp = p[idx]; m = len(p)
        c   = np.sum(1.0 / np.arange(1, m + 1))
        adj = np.minimum(1, sp * m * c / np.arange(1, m + 1))
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        out = np.empty(m); out[idx] = adj; return out

    def two_stage_bh_adj(p, alpha=0.05):
        m0  = max(n - int(np.sum(bh_adj(p) <= alpha)), 1)
        idx = np.argsort(p); sp = p[idx]
        adj = np.minimum(1, sp * m0 / np.arange(1, n + 1))
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        out = np.empty(n); out[idx] = adj; return out

    adj_bonf = np.minimum(pvals * n, 1.0)
    adj_bh   = bh_adj(pvals)
    adj_by   = by_adj(pvals)
    adj_2s   = two_stage_bh_adj(pvals)

    (
        df.select(['query_name', 'target_name', 'human_gene', 'mouse_gene',
                   'is_ortholog', 'prob_overlap', 'containment', 'jaccard'])
        .with_columns([
            pl.Series('bonferroni',   adj_bonf),
            pl.Series('bh',           adj_bh),
            pl.Series('by',           adj_by),
            pl.Series('two_stage_bh', adj_2s),
        ])
        .write_csv(f'ortholog_evaluation.hp.k{ksize}.mht.csv')
    )

    # ── Summary stats ──────────────────────────────────────────────────────────
    alpha   = 0.05
    is_orth = df['is_ortholog'].to_numpy()

    def mht_stats(adj):
        rej = adj <= alpha
        tp  = int((rej & is_orth).sum())
        tot = int(rej.sum())
        return {
            'rejected':  tot,
            'TP':        tp,
            'precision': round(tp / tot    if tot    else 0.0, 4),
            'recall':    round(tp / n_orth if n_orth else 0.0, 4),
        }

    mht_summary = {
        'bonferroni':   mht_stats(adj_bonf),
        'bh':           mht_stats(adj_bh),
        'by':           mht_stats(adj_by),
        'two_stage_bh': mht_stats(adj_2s),
    }

    metric_stats = {}
    for m in metrics:
        ov = orth_df[m].drop_nulls().to_numpy()
        nv = nort_df[m].drop_nulls().to_numpy()
        metric_stats[m] = {
            'ortholog':     {'mean': float(np.mean(ov)),   'median': float(np.median(ov)),   'n': len(ov)},
            'non_ortholog': {'mean': float(np.mean(nv)),   'median': float(np.median(nv)),   'n': len(nv)},
        }

    summary_json = {
        'ksize': ksize, 'encoding': 'hp',
        'total_hits': n_total, 'n_ortholog': n_orth, 'n_non_ortholog': n_non,
        'mht': mht_summary, 'metric_stats': metric_stats,
    }

    with open(f'ortholog_evaluation.hp.k{ksize}.summary.txt', 'w') as f:
        f.write(f'K-size: {ksize}\\nEncoding: hp\\n\\n')
        f.write(f'Total hits: {n_total:,}  |  Orthologs: {n_orth:,}  |  Non-orthologs: {n_non:,}\\n\\n')

        f.write('=== MHT Rejections (alpha=0.05) ===\\n')
        for method, s in mht_summary.items():
            f.write(f'  {method:15s}: {s["rejected"]:>10,} rejected  '
                    f'TP={s["TP"]:>8,}  prec={s["precision"]:.4f}  rec={s["recall"]:.4f}\\n')

        f.write('\\n=== Metric Means (orthologs vs non-orthologs) ===\\n')
        for m, st in metric_stats.items():
            ov, nv = st['ortholog'], st['non_ortholog']
            f.write(f'  {m}:\\n')
            f.write(f'    ortholog     mean={ov["mean"]:.4f}  median={ov["median"]:.4f}  n={ov["n"]:,}\\n')
            f.write(f'    non-ortholog mean={nv["mean"]:.4f}  median={nv["median"]:.4f}  n={nv["n"]:,}\\n')

        f.write('\\n=== JSON SUMMARY ===\\n')
        f.write(json.dumps(summary_json, indent=2) + '\\n')

    # ── Full evaluation TSV ────────────────────────────────────────────────────
    df.write_csv(f'ortholog_evaluation.hp.k{ksize}.tsv', separator='\\t')

    # ── ROC data (sampled to ≤10k points) ─────────────────────────────────────
    n_pos_t   = n_orth or 1
    n_neg_t   = n_non  or 1
    sorted_df = df.sort('prob_overlap')
    cum_tp    = sorted_df['is_ortholog'].cast(pl.Int32).cum_sum().to_numpy()
    cum_fp    = (~sorted_df['is_ortholog']).cast(pl.Int32).cum_sum().to_numpy()
    step      = max(1, n_total // 10_000)
    pl.DataFrame({
        'threshold': sorted_df['prob_overlap'].to_numpy()[::step],
        'TPR':       cum_tp[::step] / n_pos_t,
        'FPR':       cum_fp[::step] / n_neg_t,
    }).write_csv(f'ortholog_evaluation.hp.k{ksize}.roc_data.tsv', separator='\\t')
    """
}

process aggregateResults {
    publishDir params.outdir, mode: 'copy'

    input:
    path summaries

    output:
    path 'kmer_sweep_summary.tsv'
    path 'kmer_sweep_summary.json'

    script:
    """
    #!/usr/bin/env python3
    import glob, json, re

    results = []
    for f in sorted(glob.glob('*.summary.txt')):
        m = re.search(r'k(\\d+)', f)
        if not m:
            continue
        with open(f) as fh:
            content = fh.read()
        jm = re.search(r'=== JSON SUMMARY ===\\n(\\{[\\s\\S]+\\})', content)
        if jm:
            try:
                results.append(json.loads(jm.group(1)))
                continue
            except json.JSONDecodeError:
                pass
        results.append({'ksize': int(m.group(1))})

    results.sort(key=lambda x: x['ksize'])

    MHT_METHODS = ['bonferroni', 'bh', 'by', 'two_stage_bh']

    # Build TSV header
    headers = ['ksize', 'total_hits', 'n_ortholog', 'n_non_ortholog']
    for method in MHT_METHODS:
        headers += [f'{method}_rejected', f'{method}_precision', f'{method}_recall']

    with open('kmer_sweep_summary.tsv', 'w') as f:
        f.write('\\t'.join(headers) + '\\n')
        for r in results:
            row = [
                str(r.get('ksize', '')),
                str(r.get('total_hits', '')),
                str(r.get('n_ortholog', '')),
                str(r.get('n_non_ortholog', '')),
            ]
            mht = r.get('mht', {})
            for method in MHT_METHODS:
                s = mht.get(method, {})
                row += [
                    str(s.get('rejected', '')),
                    f"{s.get('precision', 0):.4f}",
                    f"{s.get('recall', 0):.4f}",
                ]
            f.write('\\t'.join(row) + '\\n')

    with open('kmer_sweep_summary.json', 'w') as f:
        json.dump({'encoding': 'hp', 'results': results}, f, indent=2)
    """
}


workflow {
    // Build kmerseek in release mode
    kmerseek_bin = buildRelease()

    // Download and parse ortholog mapping
    ortholog_file = downloadOrthologMapping()
    (ortholog_pairs, ortholog_stats) = parseOrthologMapping(ortholog_file)

    // Create k-size channel (24-40 for HP encoding)
    ksizes = Channel.of(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)

    // FASTA files are already uncompressed - use directly
    human_decompressed = channel.of(tuple('human', file(params.human_fasta)))
    mouse_decompressed = channel.of(tuple('mouse', file(params.mouse_fasta)))

    // Only index mouse (target) database - human queries are processed from FASTA
    // to avoid loading both full indices into memory simultaneously
    mouse_index_params = mouse_decompressed.combine(ksizes).map { species, fasta, ksize -> tuple(species, fasta, ksize) }

    // Index mouse database
    indexed = indexDatabase(kmerseek_bin, mouse_index_params)
    index_only = indexed[0]
    // (species, ksize, index_path)

    // Get mouse indexes by ksize
    mouse_indexes = index_only.map { species, ksize, index -> tuple(ksize, index) }

    // Combine human FASTA with mouse indexes by ksize
    // Human queries are processed from FASTA on-the-fly to save memory
    human_fasta_with_ksize = human_decompressed
        .combine(ksizes)
        .map { species, fasta, ksize -> tuple(ksize, fasta) }

    search_inputs = human_fasta_with_ksize.join(mouse_indexes)
    // (ksize, human_fasta, mouse_index)

    // Search human against mouse
    search_outputs = searchHumanVsMouse(kmerseek_bin, search_inputs)
    search_results = search_outputs[0]
    // (ksize, results_csv)

    // Evaluate ortholog detection
    eval_inputs = search_results.combine(ortholog_pairs)
    // (ksize, results_csv, ortholog_pairs)
    eval_outputs = evaluateOrthologs(eval_inputs)

    // Collect all summary files for aggregation
    // eval_outputs: [0] = (ksize, eval_tsv), [1] = summary.txt, [2] = roc_data.tsv
    summaries = eval_outputs[1].collect()

    // Aggregate results across all k-sizes
    aggregateResults(summaries)

    // Print progress
    eval_outputs[0].subscribe { ksize, eval_file ->
        println("Completed evaluation: k=${ksize} -> ${eval_file}")
    }
}
