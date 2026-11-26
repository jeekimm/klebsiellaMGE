#!/usr/bin/env python3
"""Plasmid prediction reconciliation of MOB-suite and 

`reconcile_plasmids.py` automates the table reconciliation in
`2026-11-10_plasmer_output.ipynb`. It starts from existing Plasmer and
MOB-suite calls and an existing PLSDB top-hit table. Plasmer, MOB-suite,
minimap2 and PLSDB hit selection are upstream steps.

Usage:
bash
base=/scratch/users/k2696003/klebsiella

python3 reconcile_plasmids.py \
  --base "$base" \
  --mob-tsv /path/to/mobtyper_with_contig.tsv \
  --out "$base/plasmid_reconciliation_run"


Set `--mob-tsv` to your compiled MOB-suite contig-level TSV.
The required columns are listed below.

Use a new output directory for each run. A successful run prints the
numbers of isolates with Plasmer class calls, plasmid contigs, PLSDB hits,
reconciled contigs, and contigs with replicon or relaxase evidence.

## Inputs

Paths in the table are relative to `--base`. Sample IDs are treated as
text. Keep them consistent across inputs, including suffixes such as
`v2`.

| Input | Default location or option | Format |
| --- | --- | --- |
| Cohort | `Kleb526_samples_ls.csv` or `--samples` | CSV with a `sample` column and one row per isolate |
| Plasmer output root | `plasmer/output` or `--plasmer-dir` | Each isolate has a `SAMPLE/results/` directory |
| Contig lengths | `Annotation_files/contig_lengths.tsv` or `--contig-lengths` | TSV with `sample`, `contig`, and `length` in bp. Use the GFF3-derived table from the notebook |
| PLSDB top hits | `plasmer/2026-11-16_top_plasmer_vs_PLSDB_plasmidfinder.tsv` or `--plsdb` | TSV with `sample`, `contig_id`, `query_name`, `query_length`, and `typing`. Supply the already filtered top hit per sample and contig |
| MOB-suite calls | `--mob-tsv` | Required columns are `sample`, `contig`, `length`, `rep_type(s)`, and `relaxase_type(s)`. Include `num_contigs` if available |
| Pre-extracted FASTAs | Optional `--fasta-dir` | One `SAMPLE_plasmids.fna` per isolate with a Plasmer plasmid call |

Use `--annotation-dir` to override `Annotation_files` when resolving the
default contig-length table. Use `python3 reconcile_plasmids.py --help` for the
full option list.

Each Plasmer results directory can contain these files:

| File | Expected layout |
| --- | --- |
| `SAMPLE.plasmer.predClass.tsv` | Headerless. Contig ID and class (`plasmid` or `chromosome`) |
| `SAMPLE.plasmer.predProb.tsv` | One header row. Contig ID, chromosome probability, plasmid probability |
| `SAMPLE.plasmer.predPlasmids.taxon` | Headerless. Contig ID and predicted taxon |

Missing probability or taxonomy files leave those fields blank in the
compiled output. Missing or empty class files are recorded in sample QC.
If every isolate lacks usable class calls, the run stops.

## Outputs

| File | Contents |
| --- | --- |
| `reconciled_mobsuite_plasmer_PLSDB.tsv` | Main analysis table. One row per sample and contig called as a plasmid by Plasmer (and mapped ID via PLSDB), eligible MOB-suite calls, or both |
| `plasmer_output_all.tsv` | Plasmer chromosome and plasmid calls with available probabilities, taxonomy, and GFF3-derived length |
| `plasmer_plasmid_contigs.tsv` | Plasmer plasmid calls as a two-column sample and contig list |
| `sample_qc.tsv` | One row per cohort isolate, with Plasmer file status, call counts, assembly contig count, and MOB-suite and reconciliation flags |
| `tool_intersections.tsv` | Contig counts for each observed combination of MOB-suite, Plasmer and PLSDB flags |
| `all_plasmer_plasmids.fna` | Optional combined FASTA with `sample\|contig_id` headers. Add `--fasta-dir "$base/plasmer/plasmid_sequences"` to write it |

### Main table definitions

| Column | Meaning |
| --- | --- |
| `mobsuite_prediction` | A MOB-suite row with a resolved contig ID and, where `num_contigs` is populated, a count of one |
| `plasmer_prediction` | Plasmer classified the contig as a plasmid |
| `PLSDB_mapped` | The input top-hit table supplies a query name for this Plasmer plasmid call. This flag inherits the upstream alignment filters |
| `length`, `plasmer_length`, `query_length` | Lengths from the MOB-suite input, the GFF3-derived contig table for Plasmer calls, and the supplied PLSDB hit table, respectively |
| `rep_type_presence` | A replicon type appears in MOB-suite `rep_type(s)` or in the matched PLSDB record's `typing` field |
| `relaxase_type_presence` | A relaxase type appears in MOB-suite `relaxase_type(s)` |

Blank values, dashes and common missing-value strings count as missing
marker data. PLSDB typing describes the matched reference. The replicon
flag therefore combines direct MOB-suite typing with reference-based
evidence. The presence flags record evidence available in these inputs.

## Check the results

1. Confirm that `sample_qc.tsv` has one row per expected isolate. For the
   original Bangladesh cohort, expect 526 rows.
2. Review `plasmer_status` alongside `plasmer_classified_contigs` and
   `plasmer_plasmid_contigs`. A positive classified count with zero
   plasmid calls indicates chromosome-only output. A missing or empty
   class file calls for checking the corresponding job and output files.
   The assembly contig count comes from the supplied length table.
3. Sum `n_contigs` in `tool_intersections.tsv` and compare it with the
   number of data rows in the main table. For the notebook's original
   inputs, use 1,453 Plasmer plasmid contigs, 374 PLSDB hits and 1,456
   reconciled contigs as comparison points. The script recomputes counts
   from the files provided.
4. If a join reports duplicate identifiers, inspect the listed sample
   and contig pairs in the named input. Resolve conflicting rows at
   source before rerunning, particularly competing hits in the PLSDB
   top-hit table.

The script passed a synthetic end-to-end check covering overlapping and
tool-specific calls, missing Plasmer output, marker flags and FASTA
headers. A full run against the CREATE cohort remains unverified.
"""

import argparse
from pathlib import Path
import re
from tempfile import TemporaryDirectory

import pandas as pd

KEY = ["sample", "contig_id"]
FLAGS = ["mobsuite_prediction", "plasmer_prediction", "PLSDB_mapped"]
FINAL_COLUMNS = KEY + ["query_length", "length", "plasmer_length"] + FLAGS + [
    "rep_type(s)", "typing", "relaxase_type(s)"
]


def read_table(path, columns, sep="\t"):
    frame = pd.read_csv(path, sep=sep, dtype="string", keep_default_na=False)
    require_columns(frame, columns, path)
    return frame


def require_columns(frame, columns, label):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{label}: missing columns [{', '.join(sorted(missing))}]")


def clean_keys(frame, keys):
    frame[keys] = frame[keys].astype("string").apply(lambda col: col.str.strip())
    return frame


def unique_keys(frame, keys, label):
    if (frame[keys].isna() | frame[keys].eq("")).any().any():
        raise ValueError(f"{label}: empty identifier in {keys}")
    duplicates = frame.loc[frame.duplicated(keys, keep=False), keys].drop_duplicates()
    if not duplicates.empty:
        raise ValueError(f"{label}: duplicate identifiers, e.g. {duplicates.head().to_dict('records')}")
    return frame


def present(column):
    return ~column.astype("string").fillna("").str.strip().str.lower().isin(
        ["", "-", "na", "nan", "n/a", "null", "none", "<na>"]
    )n


def read_plasmer(path, names, skiprows=0):
    if not path.exists():
        return pd.DataFrame(columns=names, dtype="string")
    try:
        frame = pd.read_csv(path, sep="\t", header=None, skiprows=skiprows,
                            dtype="string", keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=names, dtype="string")
    if len(frame.columns) != len(names):
        raise ValueError(f"{path}: expected {len(names)} columns, got {len(frame.columns)}")
    frame.columns = names
    return unique_keys(clean_keys(frame, ["contig_id"]), ["contig_id"], path)


def combine_fasta(folder, plasmids, destination):
    expected = set(plasmids[KEY].itertuples(index=False, name=None))
    observed = set()
    with destination.open("w") as output:
        for sample in sorted(plasmids["sample"].unique()):
            with (folder / f"{sample}_plasmids.fna").open() as source:
                emit = False
                for line in source:
                    if line.startswith(">"):
                        contig = line[1:].split()[0]
                        key = (sample, contig)
                        emit = key in expected
                        if emit:
                            if key in observed:
                                raise ValueError(f"Duplicate FASTA record: {sample}|{contig}")
                            observed.add(key)
                            output.write(f">{sample}|{contig}\n")
                    elif emit:
                        output.write(line.rstrip() + "\n")
    if expected != observed:
        raise ValueError(f"Missing FASTA records, e.g. {sorted(expected - observed)[:5]}")


def run(args, output):
    """Build the reconciliation tables and return a run summary."""
    annotations = args.annotation_dir or args.base / "Annotation_files"
    plasmer_dir = args.plasmer_dir or args.base / "plasmer/output"
    if not plasmer_dir.is_dir():
        raise ValueError(f"Plasmer directory does not exist: {plasmer_dir}")
    samples_path = args.samples or args.base / "Kleb526_samples_ls.csv"
    samples = clean_keys(read_table(samples_path, ["sample"], ","), ["sample"])
    cohort = unique_keys(samples, ["sample"], samples_path)["sample"]
    if cohort.empty or not cohort.map(lambda s: bool(re.fullmatch(r"[A-Za-z0-9_.-]+", s))
                                     and s not in (".", "..")).all():
        raise ValueError("Sample list must contain valid sample directory names.")

    mob_columns = ["sample", "contig", "length", "rep_type(s)", "relaxase_type(s)"]
    if args.mob_tsv:
        mob = read_table(args.mob_tsv, mob_columns)
    else:
        raise ValueError("No MOB-suite input provided.")

    mob = clean_keys(mob.rename(columns={"contig": "contig_id"}).copy(), KEY)
    mob = mob.loc[mob["sample"].isin(cohort)]
    mob_samples = set(mob["sample"])
    # Use explicit contig counts when available. The notebook's merged table
    # represents predictions without a resolved contig using missing contig IDs.
    resolved = present(mob["contig_id"])
    if "num_contigs" in mob:
        counts = pd.to_numeric(mob["num_contigs"], errors="coerce")
        resolved &= counts.eq(1) | ~present(mob["num_contigs"])
    mob = mob.loc[resolved, KEY + ["length", "rep_type(s)", "relaxase_type(s)"]].drop_duplicates()
    unique_keys(mob, KEY, "single-contig MOB-suite predictions")

    length_path = args.contig_lengths or annotations / "contig_lengths.tsv"
    lengths = read_table(length_path, ["sample", "contig", "length"]).rename(
        columns={"contig": "contig_id", "length": "plasmer_length"}
    )
    lengths = clean_keys(lengths[KEY + ["plasmer_length"]].copy(), KEY)
    lengths = lengths.loc[lengths["sample"].isin(cohort)].drop_duplicates()
    unique_keys(lengths, KEY, length_path)
    assembly_counts = lengths.groupby("sample").size()

    predictions = []
    qc = []
    for sample in cohort:
        sample_dir = plasmer_dir / sample
        prefix = sample_dir / "results" / f"{sample}.plasmer"
        class_file = Path(f"{prefix}.predClass.tsv")
        prob_file = Path(f"{prefix}.predProb.tsv")
        taxon_file = Path(f"{prefix}.predPlasmids.taxon")
        classes = read_plasmer(class_file, ["contig_id", "contig_class"])
        classes["contig_class"] = classes["contig_class"].str.strip().str.lower()
        if not classes["contig_class"].isin(["plasmid", "chromosome"]).all():
            raise ValueError(f"{class_file}: unexpected contig class")
        n_plasmids = int(classes["contig_class"].eq("plasmid").sum())
        status = "Class calls present"
        if not sample_dir.is_dir():
            status = "No output directory"
        elif not class_file.is_file():
            status = "Missing predClass.tsv"
        elif classes.empty:
            status = "Empty predClass.tsv"
        qc.append({
            "sample": sample, "plasmer_status": status,
            "predprob_present": prob_file.is_file(), "taxon_present": taxon_file.is_file(),
            "plasmer_classified_contigs": len(classes),
            "plasmer_plasmid_contigs": n_plasmids,
            "assembly_contigs": assembly_counts.get(sample, pd.NA),
            "mobsuite_report_present": sample in mob_samples,
        })
        if classes.empty:
            continue
        probs = read_plasmer(prob_file, ["contig_id", "chromosome_prob", "plasmid_prob"], 1)
        taxa = read_plasmer(taxon_file, ["contig_id", "predicted_taxon"])
        calls = classes.merge(probs, on="contig_id", how="left", validate="one_to_one")
        calls = calls.merge(taxa, on="contig_id", how="left", validate="one_to_one")
        predictions.append(calls.assign(sample=sample))

    if not predictions:
        raise ValueError("No usable Plasmer class calls found for the cohort.")
    predictions = pd.concat(predictions, ignore_index=True).merge(
        lengths, on=KEY, how="left", validate="one_to_one"
    )
    plasmids = predictions.loc[predictions["contig_class"].eq("plasmid")].copy()
    plsdb_path = args.plsdb or args.base / "plasmer/2026-11-16_top_plasmer_vs_PLSDB_plasmidfinder.tsv"
    mapping_columns = KEY + ["query_name", "query_length", "typing"]
    mapping = clean_keys(read_table(plsdb_path, mapping_columns)[mapping_columns].copy(), KEY)
    mapping = mapping.loc[mapping["sample"].isin(cohort)].drop_duplicates()
    unique_keys(mapping, KEY, plsdb_path)

    plasmids = plasmids.merge(mapping, on=KEY, how="left", validate="one_to_one")
    plasmer_final = plasmids[KEY + ["query_length", "plasmer_length", "typing"]].assign(
        plasmer_prediction=True, PLSDB_mapped=present(plasmids["query_name"])
    )
    final = plasmer_final.merge(mob.assign(mobsuite_prediction=True), on=KEY,
                                how="outer", validate="one_to_one").reindex(columns=FINAL_COLUMNS)
    for flag in FLAGS:
        final[flag] = final[flag].eq(True).fillna(False)
    for column in ["rep_type(s)", "typing", "relaxase_type(s)"]:
        final[column] = final[column].where(present(final[column]))
    final["rep_type_presence"] = final[["rep_type(s)", "typing"]].notna().any(axis=1)
    final["relaxase_type_presence"] = final["relaxase_type(s)"].notna()

    sample_qc = pd.DataFrame(qc)
    sample_qc["assembly_contigs"] = sample_qc["assembly_contigs"].astype("Int64")
    sample_qc["mobsuite_single_contig_present"] = sample_qc["sample"].isin(mob["sample"])
    sample_qc["reconciled_plasmid_present"] = sample_qc["sample"].isin(final["sample"])
    intersections = final.groupby(FLAGS).size().rename("n_contigs").reset_index()
    tables = {
        "plasmer_output_all.tsv": predictions,
        "plasmer_plasmid_contigs.tsv": plasmids[KEY],
        "reconciled_mobsuite_plasmer_PLSDB.tsv": final, # in the analysis file == combined_plasmer_results_PLSDB_mobsuite_1contig_rep_mob_types.tsv
        "sample_qc.tsv": sample_qc,
        "tool_intersections.tsv": intersections,
    }
    for name, table in tables.items():
        table.to_csv(output / name, sep="\t", index=False)
    if args.fasta_dir:
        combine_fasta(args.fasta_dir, plasmids, output / "all_plasmer_plasmids.fna")
    mapped = int(final["PLSDB_mapped"].sum())
    markers = int(final[["rep_type_presence", "relaxase_type_presence"]].any(axis=1).sum())
    return (f"Cohort: {len(cohort)} isolates. Plasmer class calls: {predictions['sample'].nunique()} isolates.\n"
            f"Plasmer plasmid contigs: {len(plasmids)}. PLSDB hits: {mapped}.\n"
            f"Reconciled plasmid contigs: {len(final)}. With replicon or relaxase evidence: {markers}.")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="reconcile_plasmids.py",
        description="Reconcile existing Plasmer, MOB-suite and top-hit PLSDB tables.",
    )
    parser.add_argument("--base", type=Path, required=True, help="Klebsiella project directory")
    parser.add_argument("--out", type=Path, required=True, help="new directory for results")
    parser.add_argument("--mob-tsv", type=Path, required=True, help="compiled MOB-suite contig-level TSV")
    for option in ("samples", "annotation-dir", "contig-lengths", "plasmer-dir", "plsdb"):
        parser.add_argument("--" + option, type=Path, help="override the project default")
    parser.add_argument("--fasta-dir", type=Path, help="optional directory with pre-extracted SAMPLE_plasmids.fna files")
    return parser


def main():
    # Stage the results so a failed analysis leaves no partial output directory.
    parser = build_parser()
    args = parser.parse_args()
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        parser.exit(1, f"Output directory already exists: {out}\n")
    try:
        with TemporaryDirectory(prefix=".reconcile-", dir=out.parent) as staging:
            staging = Path(staging)
            results = staging / "results"
            results.mkdir()
            summary = run(args, results)
            results.rename(out)
    except (OSError, ValueError, ImportError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(summary)
    print(f"Results written to {out}")


if __name__ == "__main__":
    main()
