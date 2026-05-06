# SeCAT — Sequence Consensus Amplicon Trimming

SeCAT harmonises multi-primer 16S rRNA amplicon datasets for cross-study meta-analysis. It aligns study ASV sequences to a common SILVA reference coordinate space, calculates the consensus trimming region shared across all studies, validates trim decisions against simulation-based null models, and produces a unified feature table ready for downstream community ecology analyses.

---

## How it works (brief)

1. **Clean** — Remove chloroplasts, mitochondria, empty samples; sync files.
2. **Map** — Align each study's ASVs to SILVA to find amplicon coordinates.
3. **Consensus** — Calculate the maximal overlapping 16S region across studies.
4. **Simulate** — Generate synthetic communities to build null distributions for each study's expected trimming degradation.
5. **Analyse** — Progressively trim real data and measure beta-diversity degradation at each step.
6. **Aggregate** — Compare real degradation to null distributions; apply changepoint detection and empirical p-value testing; produce per-study KEEP/EXCLUDE verdicts.
7. **Report** — Generate per-study PDF diagnostics for manual review.
8. **Trim & Merge** — Trim sequences to the consensus region, re-cluster across studies, and produce a unified feature table.
9. **Validate** — Multi-tier ecological validation of the harmonised dataset.

---

## Requirements

- [Nextflow](https://nextflow.io/) ≥ 23.04
- [Singularity](https://sylabs.io/singularity/) ≥ 3.8 **or** Docker ≥ 20
- SILVA 138 full-length aligned reference FASTA (`SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta`)

---

## Quick start

### 1. Prepare your manifest

Copy `assets/manifest_template.tsv` and fill in one row per study. Required columns:

| Column | Description |
|--------|-------------|
| `study_name` | Unique identifier for the study |
| `primer_name` | Primer pair name (e.g., `515F_806R`) |
| `asv_fasta_path` | Absolute path to ASV representative sequences FASTA |
| `asv_counts_path` | Absolute path to ASV feature/OTU table (TSV) |
| `taxonomy_path` | Absolute path to taxonomy assignments (TSV) |
| `metadata_path` | Absolute path to sample metadata (TSV) |

Optional columns: `metadata_variable`, `metadata_value` (for environment filtering).

### 2. Edit `params.yaml`

Set the two required fields:

```yaml
manifest:     "/path/to/your/secat_manifest.tsv"
reference_db: "/path/to/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta"
```

All other parameters have calibrated defaults — see the full parameter reference below.

### 3. Run

nextflow run main.nf \
    -profile slurm,singularity \
    -params-file params.yaml


**Resume after failure:**
```bash
nextflow run main.nf -profile slurm,singularity -params-file params.yaml -resume
```

---

## Site-specific HPC configuration

SeCAT uses generic resource labels so it runs on any SLURM or SGE cluster without modification. For clusters that require partition names or account billing flags, pass a custom config:

nextflow run main.nf \
    -profile slurm,singularity \
    -c conf/custom.config \
    -params-file params.yaml


For JASMIN LOTUS2, a pre-configured file is provided:
nextflow run main.nf \
    -profile slurm,singularity \
    -c conf/jasmin.config \
    -params-file params.yaml


If your data is on a non-standard filesystem (not automounted inside Singularity), add a bind mount:
```
singularity.runOptions = '--bind /path/to/your/data:/path/to/your/data'
```

---

## Pipeline workflow

| Stage | Process | Key output |
|-------|---------|------------|
| 0 | Data cleaning & QC | Filtered ASV tables, secat_manifest_clean.tsv |
| 1 | Study coordinate mapping | SILVA-aligned coordinates per study |
| 2 | Consensus region calculation | Global overlap coordinates |
| 3 | Simulation preparation | Task matrix, SILVA reference subset |
| 4 | Null model generation | Simulated degradation baselines |
| 5 | Real data trimming analysis | Study-level beta-diversity curves |
| 6 | Aggregation & verdict table | `verdict_data_all_levels.csv` |
| 7 | Report generation | Per-study PDF diagnostic reports |
| 8 | Study selection | `selected_studies_for_trim.txt` |
| 9 | Sequence trimming | `*_standardized.fasta` |
| 10 | Dataset merge | Combined feature table, taxonomy, FASTA |
| 11 | Validation | Multi-tier ecological QC report |

By default the pipeline pauses after stage 7 for manual verdict review. Set `auto_trim: true` in `params.yaml` to run stages 8–11 automatically.

### After manual review


# 1. Review output/aggregated_data/verdict_data_all_levels.csv
# 2. Edit selection_roster.txt to list studies to include
# 3. Resume trimming + merge + validation:
nextflow run main.nf \
    -profile slurm,singularity \
    -params-file params.yaml \
    -entry STANDARDIZE \
    -resume


---

## Output structure

```
output/
├── cleaned_data/               # Stage 0: filtered ASV tables and manifest
├── intermediate/               # Coordinates, alignments, simulation tasks
├── real_data_results/          # Per-study degradation curves (.rds)
├── simulation_results/         # Per-simulation degradation curves (.rds)
├── aggregated_data/
│   ├── verdict_data_all_levels.csv   # trimming verdict per study/level
│   └── master_verdict_table.csv
├── reports/
│   ├── *.pdf                         # per-study diagnostic plots
│   ├── nextflow_report.html          # Nextflow execution report
│   └── nextflow_timeline.html        # process timeline
├── standardized_datasets/
│   └── *_standardized.fasta          # trimmed ASV sequences
├── meta_analysis/
│   ├── combined_feature_table.tsv    # unified OTU/ASV table
│   ├── combined_taxonomy.tsv         # merged taxonomy
│   ├── combined_sequences.fasta      # merged FASTA
│   └── combined_metadata.tsv         # harmonised metadata
└── validation/
    └── outputs/                      # multi-tier validation results
```

---

## Full parameter reference

All parameters are set in `params.yaml`. Values shown below are the calibrated defaults optimised via systematic testing for minimal Type I error at alpha = 0.05.

### Required

| Parameter | Default | Description |
|-----------|---------|-------------|
| `manifest` | — | Path to study manifest TSV. Defines the entire study set entering the pipeline. |
| `reference_db` | — | Path to SILVA aligned FASTA. All trimming coordinates are defined relative to this reference. Changing SILVA version invalidates previously computed consensus regions. |

### Analysis mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `analysis_mode` | `"study"` | `"study"` = empirical alignment of each study's ASVs to SILVA (recommended). `"primer"` = use theoretical primer coordinates. Study mode captures study-specific biases (extraction protocol, platform effects). |

### Simulation engine

These control the synthetic community generator that builds null distributions. The null hypothesis is: "observed trimming degradation is consistent with random noise from the simulation model."

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `num_simulations` | 100 | Replicates per study. More = tighter null distributions, more stable p-values. <50 risks unreliable Type I error control. Runtime scales linearly. |
| `sim_max_silva_subset` | 10000 | Max SILVA sequences in simulation pool. Larger = more diverse null communities but more memory. Diminishing returns above ~15000. |
| `sim_use_prebuilt_subset` | true | Reuse cached SILVA subset. Set false after changing `sim_max_silva_subset`. |
| `sim_abundance_model` | `"lognormal"` | Species abundance distribution. Lognormal mimics real community dominance patterns; uniform underestimates rare-taxon trimming effects. |
| `sim_lognormal_mu` | 5 | Log-normal SAD mean (log scale). Higher = more even communities. |
| `sim_lognormal_sigma` | 2 | Log-normal SAD standard deviation. Higher = wider abundance range. |
| `sim_add_pcr_bias` | true | Simulate GC-content-dependent PCR amplification bias. Essential for realistic null distributions. |
| `sim_pcr_cycles` | 25 | Simulated PCR cycles. More cycles amplify GC bias exponentially. Match to your protocol. |
| `sim_pcr_gc_bias` | 0.65 | Strength of GC bias (0–1). Higher = stronger differential amplification. |
| `sim_pcr_optimal_gc` | 0.50 | Optimal GC fraction for PCR efficiency. Sequences near this value are preferentially amplified. |
| `sim_add_errors` | true | Simulate Illumina-like substitution and indel errors. |
| `sim_error_rate` | 0.003 | Per-base substitution rate. 0.003 = 0.3%, typical quality-filtered MiSeq. |
| `sim_indel_rate` | 0.00003 | Per-base indel rate. Extremely rare on Illumina. |
| `sim_error_position_bias` | true | Errors increase toward 3' end (Illumina quality decay). Critical for trimming-relevant null model. |
| `sim_add_chimeras` | false | Simulate chimeric sequences. Disabled by default (DADA2/Deblur already remove chimeras). |
| `sim_chimera_rate` | 0.02 | Fraction chimeric reads (only active when `sim_add_chimeras: true`). |

### Trimming

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `trim_step_mode` | `"absolute"` | `"absolute"` = fixed bp increments; `"scaled"` = increments proportional to amplicon length. Absolute gives uniform resolution; scaled normalises across primer pairs. |
| `trim_increment` | 100 | Base-pairs per trim step. Smaller = finer resolution but more steps. |
| `default_max_trim_steps` | 50 | Maximum steps per study. Caps the search space. |
| `consensus_buffer_steps` | 20 | Extra steps beyond consensus boundary. Ensures degradation curves have sufficient post-change data. |
| `max_absolute_trim_steps` | 2000 | Hard safety ceiling. Should never be reached normally. |

### Quality filters

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `min_relative_abundance` | 0 | Minimum relative abundance to retain a taxon. 0 = keep all. Higher values reduce noise but may mask rare-taxon impacts. |
| `min_taxa_for_bray` | 3 | Minimum taxa per sample for Bray-Curtis computation. Prevents degenerate distances from near-empty samples. |
| `min_median_read_depth` | 100 | Minimum median reads per study. Filters severely undersequenced datasets. |

### Changepoint detection (PELT)

Identifies the trim step where beta-diversity shows a statistically significant shift.

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `changepoint_penalty_method` | `"MANUAL"` | PELT penalty type. `"MANUAL"` = scaled by multiplier below. Other options: `"BIC"`, `"AIC"`. |
| `changepoint_penalty_multiplier` | 1 | Multiplier for MANUAL penalty (penalty = multiplier × variance). Higher = fewer changepoints detected (more conservative). Lower = more sensitive but higher false positive rate. `params.yaml` may override this to 50 for production runs. |

### Null model comparison

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `null_model_p_threshold` | 0.05 | P-value threshold for null rejection. Lower = more conservative (fewer false positives). |
| `null_model_min_consecutive` | 3 | Required consecutive significant steps to trigger a verdict. Reduces single-point noise. Higher = fewer false alarms but delayed detection. |
| `null_model_min_trim_bp` | 5 | Minimum bp trimmed before testing begins. Creates a "dead zone" where trivial trims aren't evaluated. |
| `distance_cutoff_threshold` | 0.15 | Absolute Bray-Curtis threshold (secondary method). Studies exceeding this are flagged regardless of null model results. Lower = stricter. |
| `distance_cutoff_min_trim_bp` | 5 | Minimum bp before distance cutoff is evaluated. |

### Consensus region

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `consensus_optimization_threshold` | 0.20 | Maximum tolerable degradation for including a study in the consensus. Higher = more inclusive but accepts more distortion. |
| `min_consensus_studies` | 3 | Minimum studies contributing to the consensus region. |
| `min_consensus_coverage` | 0.50 | Minimum fraction of input studies representable in the consensus. |

### Reference alignment

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `reference_alignment_mode` | `"subset"` | `"subset"` = random SILVA subset for speed; `"full"` = entire SILVA DB (requires 64+ GB RAM). |
| `reference_subset_size` | 5000 | Number of SILVA sequences in subset. Larger = more accurate alignment, diminishing returns above ~10000. |

### Merge

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `merge_method` | `"advanced"` | `"advanced"` = taxonomy-aware with conflict resolution; `"simple"` = direct concatenation. Advanced handles ASV identity conflicts across studies. |
| `harmonize_metadata` | true | Standardise metadata column names using synonym mapping (e.g., "lat" → "Latitude"). |
| `selection_mode` | `"roster"` | `"auto"` = include all KEEP studies; `"roster"` = read from a curated file. Manual review recommended for new datasets. |

### Pipeline control

| Parameter | Default | Description & data impact |
|-----------|---------|---------------------------|
| `auto_trim` | false | If true, skip manual review and run trim/merge immediately after verdicts. Recommended only for previously vetted datasets. |
| `run_validation` | true | Run multi-tier ecological validation after merge. |
| `validation_levels` | `"ASV,Genus,Family"` | Taxonomic levels for validation. More levels = more thorough but slower. |

---

## Resource requirements

| Process | Typical memory | CPUs | Time |
|---------|---------------|------|------|
| Data cleaning | 4 GB | 1 | < 1 h |
| Study mapping | 16–64 GB | 1–4 | 2–8 h |
| Simulation (per replicate) | 8 GB | 1 | 1–4 h |
| Real data analysis | 16 GB | 1 | 2–4 h |
| Aggregation | 16–30 GB | 1 | 2–8 h |
| Merge | 16–48 GB | 1–16 | 4–12 h |
| Validation | 8–16 GB | 1 | 1–4 h |

Total runtime depends on the number of studies and simulations. A typical 15-study run with 100 simulations completes in 12–24 hours on a SLURM cluster.

---

## Troubleshooting

**Pipeline fails at STUDY_MAPPING with memory errors:** Increase `reference_subset_size` cautiously, or switch to `reference_alignment_mode: "full"` with a high-memory node. Alternatively, set `use_all_asvs: false` and reduce `asv_sample_size`.

**Too many studies flagged as EXCLUDE:** Consider increasing `changepoint_penalty_multiplier` (makes changepoint detection more conservative), or increasing `distance_cutoff_threshold` (tolerates more compositional change).

**Too few studies flagged (permissive verdicts):** Decrease `changepoint_penalty_multiplier` or lower `null_model_p_threshold`.

**Simulations take too long:** Reduce `num_simulations` to 50 for exploratory runs. Increase back to 100+ for final analysis.

**Container pull fails:** Ensure network access to ghcr.io. Alternatively, pre-pull: `singularity pull docker://ghcr.io/derbydt/secat:latest`

---

## Citation

Terrey D. et al. (in preparation). SeCAT: a bioinformatics pipeline for cross-study harmonisation of 16S rRNA amplicon datasets.

---

## Licence

MIT
