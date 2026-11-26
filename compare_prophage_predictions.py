#!/usr/bin/env python3
"""Compare PHASTEST, geNomad and PhiSpy prophage predictions.

Inputs are preprocessed TSV tables, not raw tool output (i.e., intact phages selected from result files of PHASTEST, concatenated '*_prophage.tsv' for phispy, and '*_provirus.tsv' for geNomad).
Coordinates must refer to the same assembly, use identical contig identifiers, and be 1-based
inclusive. PHASTEST reported lengths are interpreted as kb; other lengths as bp.
Requires pandas and networkx; for the default UpSet PNG also requires
matplotlib and upsetplot (pip install pandas networkx matplotlib upsetplot).
Consensus loci use transitive connected components and enclosing spans.

Usage:
 python compare_prophage_predictions.py --phastest phastest.tsv \
   --genomad genomad.tsv --phispy phispy.tsv --output-dir results
"""
import argparse
from itertools import combinations
from pathlib import Path
import sys
import pandas as pd
import networkx as nx

SCHEMAS = {
    'PHASTEST': {'Sample': 'sample', 'Contig': 'contig', 'start': 'start', 'end': 'end', 'length': 'length'},
    'geNomad': {'sample': 'sample', 'source_seq': 'contig', 'start': 'start', 'end': 'end', 'length': 'length'},
    'PhiSpy': {'sample': 'sample', 'Contig': 'contig', 'Start': 'start', 'Stop': 'end', 'Length': 'length'},
}
MATCH_COLUMNS = ['sample', 'contig', 'tool1', 'tool2', 'start1', 'end1', 'start2', 'end2', 'reciprocal_overlap1', 'reciprocal_overlap2', 'iou']


def load_predictions(path, tool):
    table = pd.read_csv(path, sep='\t', dtype=str)
    mapping = SCHEMAS[tool]
    missing = set(mapping).difference(table.columns)
    if missing:
        raise ValueError(f'{tool}: missing required columns: {sorted(missing)}')
    df = table[list(mapping)].rename(columns=mapping).copy()
    for field in ['sample', 'contig']:
        df[field] = df[field].astype('string').str.strip()
        if df[field].isna().any() or df[field].eq('').any():
            raise ValueError(f'{tool}: empty {field} identifiers found')
    for field in ['start', 'end']:
        values = pd.to_numeric(df[field], errors='coerce')
        if values.isna().any() or (~(values % 1 == 0)).any():
            raise ValueError(f'{tool}: invalid or non-integer {field} coordinates')
        df[field] = values.astype('int64')
    if (df['start'] < 1).any() or (df['end'] < df['start']).any():
        raise ValueError(f'{tool}: coordinates must satisfy 1 <= start <= end')
    lengths = df['length'].astype('string').str.strip()
    if tool == 'PHASTEST':
        lengths = lengths.str.replace(r'\s*kb\s*$', '', case=False, regex=True)
    df['length'] = pd.to_numeric(lengths, errors='coerce')
    if df['length'].isna().any() or (df['length'] <= 0).any():
        raise ValueError(f'{tool}: invalid or missing lengths')
    if tool == 'PHASTEST':
        df['length'] *= 1000
    df['tool'] = tool
    return df[['sample', 'contig', 'tool', 'start', 'end', 'length']]


def reciprocal_overlap(a, b):
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    len_a = a.end - a.start + 1
    len_b = b.end - b.start + 1
    return overlap / len_a, overlap / len_b, overlap


def find_matches(df, threshold):
    matches = []
    for (sample, contig), group in df.groupby(['sample', 'contig'], sort=True):
        rows = list(group.itertuples(index=False))
        for a, b in combinations(rows, 2):
            if a.tool == b.tool:
                continue
            ro_a, ro_b, intersection = reciprocal_overlap(a, b)
            if ro_a >= threshold and ro_b >= threshold:
                # IoU uses the union of covered bases, not the enclosing span.
                union = (a.end - a.start + 1) + (b.end - b.start + 1) - intersection
                matches.append((sample, contig, a.tool, b.tool, a.start, a.end,
                                b.start, b.end, ro_a, ro_b, intersection / union))
    return pd.DataFrame(matches, columns=MATCH_COLUMNS)



def build_consensus_loci(df, threshold):
    """Group predictions by connected components of qualifying overlap edges.

    A component may include indirectly linked predictions. Its start/end
    describe the enclosing span, not an intersection or precise boundary.
    Isolated predictions become single-tool loci.
    """
    graph = nx.Graph()
    graph.add_nodes_from(df.index)
    # Match records contain genomic coordinates but not row indices. Rebuild
    # qualifying edges using the same interval test and row indices.
    for (_, _), group in df.groupby(['sample', 'contig'], sort=True):
        records = list(group.iterrows())
        for (i, a), (j, b) in combinations(records, 2):
            if a['tool'] == b['tool']:
                continue
            overlap = max(0, min(a['end'], b['end']) - max(a['start'], b['start']) + 1)
            if overlap == 0:
                continue
            if (overlap / (a['end'] - a['start'] + 1) >= threshold
                    and overlap / (b['end'] - b['start'] + 1) >= threshold):
                graph.add_edge(i, j)

    rows = []
    components = sorted(nx.connected_components(graph), key=lambda c: min(c))
    for locus_number, component in enumerate(components, start=1):
        subset = df.loc[sorted(component)]
        tools = sorted(subset['tool'].unique())
        rows.append({
            'locus_id': f'locus_{locus_number:06d}',
            'sample': subset['sample'].iloc[0],
            'contig': subset['contig'].iloc[0],
            'start': int(subset['start'].min()),
            'end': int(subset['end'].max()),
            'tools': ','.join(tools),
            'tool_combination': ' + '.join(tools),
            'n_tools': len(tools),
            'n_predictions': len(subset),
            'enclosing_span_bp': int(subset['end'].max() - subset['start'].min() + 1),
        })
    columns = ['locus_id', 'sample', 'contig', 'start', 'end', 'tools',
               'tool_combination', 'n_tools', 'n_predictions', 'enclosing_span_bp']
    return pd.DataFrame(rows, columns=columns)


def save_upset_plot(loci, path):
    """Save an UpSet plot of distinct consensus loci by tool membership."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_memberships

    if loci.empty:
        print('No loci found; skipping UpSet plot.')
        return
    memberships = [item.split(',') for item in loci['tools']]
    upset_data = from_memberships(memberships, data=[1] * len(memberships))
    fig = plt.figure(figsize=(11, 6))
    upset = UpSet(upset_data, subset_size='sum', show_counts=True,
                  sort_by='cardinality', sort_categories_by=None)
    upset.plot(fig=fig)
    fig.suptitle('Prophage detection agreement between PHASTEST, geNomad and PhiSpy',
                 fontsize=14)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def summarize(df):
    counts = df.groupby(['tool', 'sample']).size().rename('prophage_count').reset_index()
    stats = (counts.groupby('tool')['prophage_count']
             .agg(total_isolates_with_prophages='count', mean='mean', median='median',
                  std='std', min='min', max='max').T)
    stats.index.name = 'statistic'
    return stats.reset_index(), counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--phastest', type=Path, required=True, help='Preprocessed PHASTEST intact-prophage TSV')
    parser.add_argument('--genomad', type=Path, required=True, help='Preprocessed geNomad provirus TSV')
    parser.add_argument('--phispy', type=Path, required=True, help='Preprocessed PhiSpy prophage TSV')
    parser.add_argument('--output-dir', type=Path, required=True, help='Output directory')
    parser.add_argument('--no-upset', action='store_true', help='Skip generating the UpSet PNG')
    parser.add_argument('--threshold', type=float, default=0.5, help='Minimum reciprocal overlap for both predictions (default: 0.5)')
    args = parser.parse_args()
    if not 0 <= args.threshold <= 1:
        parser.error('--threshold must be between 0 and 1')
    combined = pd.concat([load_predictions(args.phastest, 'PHASTEST'),
                          load_predictions(args.genomad, 'geNomad'),
                          load_predictions(args.phispy, 'PhiSpy')], ignore_index=True)
    stats, counts = summarize(combined)
    matches = find_matches(combined, args.threshold)
    loci = build_consensus_loci(combined, args.threshold)
    combinations_summary = (loci['tool_combination'].value_counts()
                            .rename_axis('tool_combination').reset_index(name='n_loci'))
    support_summary = (loci['n_tools'].value_counts().sort_index()
                       .rename_axis('n_tools').reset_index(name='n_loci'))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in [('combined_prophage_predictions.tsv', combined),
                        ('prophage_summary_statistics.tsv', stats),
                        ('prophage_counts_per_sample.tsv', counts),
                        ('qualifying_pairwise_overlaps.tsv', matches),
                        ('prophage_consensus_loci.tsv', loci),
                        ('prophage_tool_combination_counts.tsv', combinations_summary),
                        ('prophage_tool_support_counts.tsv', support_summary)]:
        frame.to_csv(args.output_dir / name, sep='\t', index=False)
    if not args.no_upset:
        save_upset_plot(loci, args.output_dir / 'prophage_upset_plot.png')
    print(f'Predictions: {len(combined)}; qualifying pairwise matches: {len(matches)}; consensus loci: {len(loci)}')
    print(f'Output directory: {args.output_dir.resolve()}')
    print(stats.to_string(index=False))
    print('\nLoci by tool combination:')
    print(combinations_summary.to_string(index=False))
    print('\nLoci by number of supporting tools:')
    print(support_summary.to_string(index=False))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
