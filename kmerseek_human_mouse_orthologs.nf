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

// K-mer size range for HP encoding (15-30)
params.k_min = 15
params.k_max = 30

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

    script:
    """
    #!/usr/bin/env python3
    import csv
    import re
    import json
    from collections import defaultdict

    # Load ortholog pairs (ground truth)
    ortholog_pairs_dict = defaultdict(set)
    with open('${ortholog_pairs}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            human_gene = row['human_gene'].upper()
            mouse_gene = row['mouse_gene']
            ortholog_pairs_dict[human_gene].add(mouse_gene)

    # Extract gene symbol from GENCODE FASTA header
    # Format: ENSP...|ENST...|ENSG...|...|...|SYMBOL-201|SYMBOL|length
    def extract_gene_symbol(header):
        parts = header.split('|')
        if len(parts) >= 7:
            return parts[6]  # Gene symbol is 7th field (0-indexed: 6)
        return None

    # Process kmerseek results
    # Expected columns: query_name, target_name, score, etc.
    results = []
    with open('${results_csv}', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            query_header = row['query_name']
            target_header = row['target_name']

            human_gene = extract_gene_symbol(query_header)
            mouse_gene = extract_gene_symbol(target_header)

            if human_gene and mouse_gene:
                # Check if human gene is uppercase (as expected)
                human_gene_upper = human_gene.upper()

                # Check if this is a true ortholog pair
                is_ortholog = mouse_gene in ortholog_pairs_dict.get(human_gene_upper, set())

                # Get score for ranking (try multiple possible score columns)
                score = 0.0
                for score_col in ['average_kmer_rarity', 'score', 'jaccard', 'containment']:
                    if score_col in row and row[score_col]:
                        try:
                            score = float(row[score_col])
                            break
                        except ValueError:
                            continue

                results.append({
                    'query_header': query_header,
                    'target_header': target_header,
                    'human_gene': human_gene,
                    'mouse_gene': mouse_gene,
                    'is_ortholog': is_ortholog,
                    'score': score,
                    **{k: v for k, v in row.items() if k not in ['query_name', 'target_name']}
                })

    # Write detailed evaluation results
    with open('ortholog_evaluation.hp.k${ksize}.tsv', 'w') as f:
        if results:
            fieldnames = ['human_gene', 'mouse_gene', 'is_ortholog', 'score'] + \\
                        [k for k in results[0].keys() if k not in ['human_gene', 'mouse_gene', 'is_ortholog', 'score']]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\\t')
            writer.writeheader()
            for row in results:
                writer.writerow(row)

    # ============================================
    # ML METRICS CALCULATION
    # ============================================

    # Hit-level metrics
    total_hits = len(results)
    TP_hits = sum(1 for r in results if r['is_ortholog'])  # True positives at hit level
    FP_hits = total_hits - TP_hits  # False positives at hit level

    # Precision at hit level
    precision_hits = TP_hits / total_hits if total_hits > 0 else 0

    # Gene-level analysis
    # Get all unique human genes that were queried
    human_genes_searched = set(r['human_gene'].upper() for r in results)
    human_genes_in_ground_truth = set(ortholog_pairs_dict.keys())

    # For recall: count unique ortholog pairs found vs total possible
    found_pairs = set()
    for r in results:
        if r['is_ortholog']:
            found_pairs.add((r['human_gene'].upper(), r['mouse_gene']))

    # Total ortholog pairs that could have been found (human gene was searched)
    searchable_pairs = set()
    for human_gene in human_genes_searched:
        for mouse_gene in ortholog_pairs_dict.get(human_gene, set()):
            searchable_pairs.add((human_gene, mouse_gene))

    TP_pairs = len(found_pairs)
    FN_pairs = len(searchable_pairs) - TP_pairs

    # Recall (sensitivity) at pair level
    sensitivity = TP_pairs / len(searchable_pairs) if searchable_pairs else 0
    recall = sensitivity  # Alias

    # For specificity, we need to consider negative pairs
    # Negative pairs = pairs that are NOT orthologs but were returned as hits
    FP_pairs = sum(1 for r in results if not r['is_ortholog'])

    # Gene-level recall: fraction of human genes that found at least one correct ortholog
    human_genes_with_orthologs_found = set(p[0] for p in found_pairs)
    searchable_human_genes = human_genes_searched & human_genes_in_ground_truth
    gene_level_recall = len(human_genes_with_orthologs_found) / len(searchable_human_genes) if searchable_human_genes else 0

    # F1 score
    f1 = 2 * precision_hits * recall / (precision_hits + recall) if (precision_hits + recall) > 0 else 0

    # ============================================
    # BEST HIT ANALYSIS
    # ============================================
    best_hits = {}
    for result in results:
        query = result['query_header']
        score = result['score']
        if query not in best_hits or score > best_hits[query]['score']:
            best_hits[query] = {
                'human_gene': result['human_gene'],
                'mouse_gene': result['mouse_gene'],
                'is_ortholog': result['is_ortholog'],
                'score': score
            }

    best_hit_TP = sum(1 for h in best_hits.values() if h['is_ortholog'])
    best_hit_FP = len(best_hits) - best_hit_TP
    best_hit_precision = best_hit_TP / len(best_hits) if best_hits else 0

    # ============================================
    # ROC CURVE AND AUC CALCULATION
    # ============================================
    # Sort results by score (descending) for ROC curve
    sorted_results = sorted(results, key=lambda x: -x['score'])

    # Calculate ROC curve points
    roc_points = []
    cumulative_TP = 0
    cumulative_FP = 0
    total_positives = TP_hits
    total_negatives = FP_hits

    prev_score = None
    for i, r in enumerate(sorted_results):
        if r['is_ortholog']:
            cumulative_TP += 1
        else:
            cumulative_FP += 1

        # Calculate TPR (sensitivity) and FPR (1 - specificity)
        TPR = cumulative_TP / total_positives if total_positives > 0 else 0
        FPR = cumulative_FP / total_negatives if total_negatives > 0 else 0

        # Add point when score changes or at the end
        if prev_score is None or r['score'] != prev_score or i == len(sorted_results) - 1:
            roc_points.append({
                'threshold': r['score'],
                'TPR': TPR,
                'FPR': FPR,
                'cumulative_TP': cumulative_TP,
                'cumulative_FP': cumulative_FP
            })
        prev_score = r['score']

    # Calculate AUC using trapezoidal rule
    auc = 0.0
    for i in range(1, len(roc_points)):
        # Trapezoid area = (x2 - x1) * (y1 + y2) / 2
        dx = roc_points[i]['FPR'] - roc_points[i-1]['FPR']
        avg_y = (roc_points[i]['TPR'] + roc_points[i-1]['TPR']) / 2
        auc += dx * avg_y

    # Write ROC data
    with open('ortholog_evaluation.hp.k${ksize}.roc_data.tsv', 'w') as f:
        f.write('threshold\\tTPR\\tFPR\\tcumulative_TP\\tcumulative_FP\\n')
        for point in roc_points:
            f.write(f"{point['threshold']:.6f}\\t{point['TPR']:.6f}\\t{point['FPR']:.6f}\\t{point['cumulative_TP']}\\t{point['cumulative_FP']}\\n")

    # ============================================
    # PRECISION-RECALL CURVE
    # ============================================
    pr_points = []
    cumulative_TP = 0
    for i, r in enumerate(sorted_results):
        if r['is_ortholog']:
            cumulative_TP += 1
        total_predicted = i + 1
        prec = cumulative_TP / total_predicted
        rec = cumulative_TP / total_positives if total_positives > 0 else 0
        pr_points.append({'precision': prec, 'recall': rec, 'threshold': r['score']})

    # Calculate Average Precision (AP)
    ap = 0.0
    prev_recall = 0
    for point in pr_points:
        if point['recall'] > prev_recall:
            ap += point['precision'] * (point['recall'] - prev_recall)
            prev_recall = point['recall']

    # ============================================
    # ADDITIONAL METRICS AT VARIOUS THRESHOLDS
    # ============================================
    # Find metrics at specific recall levels (0.9, 0.8, 0.7, etc.)
    recall_thresholds = [0.9, 0.8, 0.7, 0.6, 0.5]
    precision_at_recall = {}
    for target_recall in recall_thresholds:
        for point in pr_points:
            if point['recall'] >= target_recall:
                precision_at_recall[target_recall] = point['precision']
                break
        else:
            precision_at_recall[target_recall] = 0.0

    # ============================================
    # TOP-K ACCURACY
    # ============================================
    # For each query, check if ortholog is in top-K results
    query_results = defaultdict(list)
    for r in results:
        query_results[r['query_header']].append(r)

    top_k_values = [1, 3, 5, 10, 20]
    top_k_accuracy = {}
    for k in top_k_values:
        correct = 0
        total_queries_with_orthologs = 0
        for query, hits in query_results.items():
            # Sort hits by score
            sorted_hits = sorted(hits, key=lambda x: -x['score'])
            # Check if query's gene has known orthologs
            human_gene = sorted_hits[0]['human_gene'].upper()
            if human_gene in ortholog_pairs_dict:
                total_queries_with_orthologs += 1
                # Check if any of top-K are orthologs
                top_k_hits = sorted_hits[:k]
                if any(h['is_ortholog'] for h in top_k_hits):
                    correct += 1
        top_k_accuracy[k] = correct / total_queries_with_orthologs if total_queries_with_orthologs > 0 else 0

    # ============================================
    # MEAN RECIPROCAL RANK (MRR)
    # ============================================
    reciprocal_ranks = []
    for query, hits in query_results.items():
        sorted_hits = sorted(hits, key=lambda x: -x['score'])
        human_gene = sorted_hits[0]['human_gene'].upper()
        if human_gene in ortholog_pairs_dict:
            # Find rank of first correct ortholog
            for rank, hit in enumerate(sorted_hits, start=1):
                if hit['is_ortholog']:
                    reciprocal_ranks.append(1.0 / rank)
                    break
            else:
                reciprocal_ranks.append(0.0)  # No ortholog found

    mrr = sum(reciprocal_ranks) / len(reciprocal_ranks) if reciprocal_ranks else 0

    # ============================================
    # WRITE SUMMARY
    # ============================================
    with open('ortholog_evaluation.hp.k${ksize}.summary.txt', 'w') as f:
        f.write(f'K-size: ${ksize}\\n')
        f.write(f'Encoding: hp\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('HIT-LEVEL METRICS\\n')
        f.write('='*50 + '\\n')
        f.write(f'Total search hits: {total_hits}\\n')
        f.write(f'True positive hits (orthologs): {TP_hits}\\n')
        f.write(f'False positive hits (non-orthologs): {FP_hits}\\n')
        f.write(f'Precision: {precision_hits:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('PAIR-LEVEL METRICS\\n')
        f.write('='*50 + '\\n')
        f.write(f'Searchable ortholog pairs: {len(searchable_pairs)}\\n')
        f.write(f'True positives (pairs found): {TP_pairs}\\n')
        f.write(f'False negatives (pairs missed): {FN_pairs}\\n')
        f.write(f'Sensitivity (Recall): {sensitivity:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('GENE-LEVEL METRICS\\n')
        f.write('='*50 + '\\n')
        f.write(f'Human genes searched: {len(human_genes_searched)}\\n')
        f.write(f'Human genes with known orthologs: {len(human_genes_in_ground_truth)}\\n')
        f.write(f'Human genes searched with known orthologs: {len(searchable_human_genes)}\\n')
        f.write(f'Human genes with correct ortholog found: {len(human_genes_with_orthologs_found)}\\n')
        f.write(f'Gene-level recall: {gene_level_recall:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('COMPOSITE METRICS\\n')
        f.write('='*50 + '\\n')
        f.write(f'F1 Score: {f1:.4f}\\n')
        f.write(f'AUC-ROC: {auc:.4f}\\n')
        f.write(f'Average Precision (AP): {ap:.4f}\\n')
        f.write(f'Mean Reciprocal Rank (MRR): {mrr:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('BEST HIT ANALYSIS\\n')
        f.write('='*50 + '\\n')
        f.write(f'Total unique queries: {len(best_hits)}\\n')
        f.write(f'Best hits that are orthologs (TP): {best_hit_TP}\\n')
        f.write(f'Best hits that are not orthologs (FP): {best_hit_FP}\\n')
        f.write(f'Best hit precision: {best_hit_precision:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('TOP-K ACCURACY\\n')
        f.write('='*50 + '\\n')
        for k in top_k_values:
            f.write(f'Top-{k} accuracy: {top_k_accuracy[k]:.4f}\\n')
        f.write(f'\\n')

        f.write('='*50 + '\\n')
        f.write('PRECISION AT RECALL LEVELS\\n')
        f.write('='*50 + '\\n')
        for r, p in sorted(precision_at_recall.items(), reverse=True):
            f.write(f'Precision @ Recall={r:.1f}: {p:.4f}\\n')
        f.write(f'\\n')

        # JSON summary for easy parsing
        f.write('='*50 + '\\n')
        f.write('JSON SUMMARY\\n')
        f.write('='*50 + '\\n')
        summary_json = {
            'ksize': ${ksize},
            'encoding': 'hp',
            'total_hits': total_hits,
            'TP_hits': TP_hits,
            'FP_hits': FP_hits,
            'precision': round(precision_hits, 4),
            'sensitivity': round(sensitivity, 4),
            'recall': round(recall, 4),
            'gene_level_recall': round(gene_level_recall, 4),
            'f1': round(f1, 4),
            'auc_roc': round(auc, 4),
            'average_precision': round(ap, 4),
            'mrr': round(mrr, 4),
            'best_hit_precision': round(best_hit_precision, 4),
            'top_k_accuracy': {f'top_{k}': round(v, 4) for k, v in top_k_accuracy.items()},
            'precision_at_recall': {f'recall_{r}': round(p, 4) for r, p in precision_at_recall.items()}
        }
        f.write(json.dumps(summary_json, indent=2))
        f.write('\\n')
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
    import re
    import os
    import glob
    import json

    results = []

    for summary_file in glob.glob('*.summary.txt'):
        ksize_match = re.search(r'k(\\d+)', summary_file)
        if not ksize_match:
            continue
        ksize = int(ksize_match.group(1))

        with open(summary_file) as f:
            content = f.read()

        # Try to parse JSON summary from the file
        json_match = re.search(r'JSON SUMMARY\\n=+\\n(\\{[\\s\\S]+\\})', content)
        if json_match:
            try:
                metrics = json.loads(json_match.group(1))
                results.append(metrics)
                continue
            except json.JSONDecodeError:
                pass

        # Fallback: parse metrics manually
        def extract_float(pattern, default=0.0):
            match = re.search(pattern, content)
            return float(match.group(1)) if match else default

        def extract_int(pattern, default=0):
            match = re.search(pattern, content)
            return int(match.group(1)) if match else default

        results.append({
            'ksize': ksize,
            'precision': extract_float(r'Precision: ([\\d.]+)'),
            'sensitivity': extract_float(r'Sensitivity \\(Recall\\): ([\\d.]+)'),
            'recall': extract_float(r'Sensitivity \\(Recall\\): ([\\d.]+)'),
            'gene_level_recall': extract_float(r'Gene-level recall: ([\\d.]+)'),
            'f1': extract_float(r'F1 Score: ([\\d.]+)'),
            'auc_roc': extract_float(r'AUC-ROC: ([\\d.]+)'),
            'average_precision': extract_float(r'Average Precision \\(AP\\): ([\\d.]+)'),
            'mrr': extract_float(r'Mean Reciprocal Rank \\(MRR\\): ([\\d.]+)'),
            'best_hit_precision': extract_float(r'Best hit precision: ([\\d.]+)'),
            'total_hits': extract_int(r'Total search hits: (\\d+)'),
            'TP_hits': extract_int(r'True positive hits \\(orthologs\\): (\\d+)'),
            'FP_hits': extract_int(r'False positive hits \\(non-orthologs\\): (\\d+)'),
            'top_k_accuracy': {
                'top_1': extract_float(r'Top-1 accuracy: ([\\d.]+)'),
                'top_3': extract_float(r'Top-3 accuracy: ([\\d.]+)'),
                'top_5': extract_float(r'Top-5 accuracy: ([\\d.]+)'),
                'top_10': extract_float(r'Top-10 accuracy: ([\\d.]+)'),
                'top_20': extract_float(r'Top-20 accuracy: ([\\d.]+)')
            }
        })

    # Sort by ksize
    results.sort(key=lambda x: x['ksize'])

    # Write comprehensive TSV summary table
    with open('kmer_sweep_summary.tsv', 'w') as f:
        headers = [
            'ksize', 'precision', 'sensitivity', 'gene_level_recall', 'f1',
            'auc_roc', 'average_precision', 'mrr', 'best_hit_precision',
            'top_1_accuracy', 'top_3_accuracy', 'top_5_accuracy', 'top_10_accuracy', 'top_20_accuracy',
            'total_hits', 'TP_hits', 'FP_hits'
        ]
        f.write('\\t'.join(headers) + '\\n')

        for r in results:
            top_k = r.get('top_k_accuracy', {})
            row = [
                str(r.get('ksize', '')),
                f"{r.get('precision', 0):.4f}",
                f"{r.get('sensitivity', 0):.4f}",
                f"{r.get('gene_level_recall', 0):.4f}",
                f"{r.get('f1', 0):.4f}",
                f"{r.get('auc_roc', 0):.4f}",
                f"{r.get('average_precision', 0):.4f}",
                f"{r.get('mrr', 0):.4f}",
                f"{r.get('best_hit_precision', 0):.4f}",
                f"{top_k.get('top_1', 0):.4f}",
                f"{top_k.get('top_3', 0):.4f}",
                f"{top_k.get('top_5', 0):.4f}",
                f"{top_k.get('top_10', 0):.4f}",
                f"{top_k.get('top_20', 0):.4f}",
                str(r.get('total_hits', 0)),
                str(r.get('TP_hits', 0)),
                str(r.get('FP_hits', 0))
            ]
            f.write('\\t'.join(row) + '\\n')

    # Write JSON for programmatic access
    with open('kmer_sweep_summary.json', 'w') as f:
        json.dump({
            'encoding': 'hp',
            'k_range': [min(r['ksize'] for r in results), max(r['ksize'] for r in results)],
            'results': results
        }, f, indent=2)
    """
}

workflow {
    // Build kmerseek in release mode
    kmerseek_bin = buildRelease()

    // Download and parse ortholog mapping
    ortholog_file = downloadOrthologMapping()
    (ortholog_pairs, ortholog_stats) = parseOrthologMapping(ortholog_file)

    // Create k-size channel (15-30 for HP encoding)
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
