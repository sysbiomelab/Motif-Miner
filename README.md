# Motif-Miner

Nextflow pipeline for finding conserved regulatory motifs upstream of genes in families of biosynthetic gene clusters (BGCs).

Given BGC families as antiSMASH GenBank files, the pipeline groups the proteins of each family into orthologous clusters (COGs), extracts the upstream intergenic region of every gene in each COG, removes redundant sequences, runs MEME on each COG, and sorts the resulting motifs by the antiSMASH gene kind of the genes they sit in front of. It was written for the analysis in the paper cited below.

## How it works

Each family is processed independently:

1. GenBank files are converted to protein FASTA and a CDS table is built (locus tag, coordinates, product, antiSMASH `gene_kind`).
2. ProteinOrtho groups the proteins into COGs. One FASTA per COG is written, numbered in the order of the ProteinOrtho table.
3. For every gene in every COG, the region between the previous same-strand CDS (or the contig edge) and the gene start is extracted, plus the first `inside_gene_len` bp of the gene. Regions are capped at `max_upstream_len`, reverse-complemented on the minus strand, and dropped if shorter than `min_upstream_len`.
4. CD-HIT-EST removes near-identical upstream sequences within each COG.
5. The percentage of regulatory, biosynthetic and biosynthetic-additional genes in each COG is computed from the CDS table.
6. MEME is run on every COG with more than four remaining sequences, six times with different settings (see below). Families with no MEME output are excluded and listed in `families_excluded_no_meme_results.txt`.
7. MEME XML is parsed into CSV tables and MEME-format PWM files.
8. Each COG's results are copied into `regulatory/`, `biosynthetic/`, `biosynthetic_additional/` and/or `unknown/` according to its gene kind percentages and `org_threshold`.
9. A per-family table joins every motif site to its gene annotation.

## Requirements

- Nextflow 21.10 or later, Java 11 or later
- Singularity/Apptainer, or Docker with `-profile docker`, for ProteinOrtho 6.3.5, CD-HIT 4.8.1 and MEME 5.5.8 (pulled as version-pinned containers)
- Conda or Mamba for the Python steps (Python 3.9, Biopython 1.84, pandas 2.2; the environment is built from `conf/conda.yml` on first run)

## Input

The input directory holds one subdirectory per family, named `FAM*`, each containing the antiSMASH GenBank files (`.gbk` or `.gb`, one region per file) of that family's BGCs:

```
bgc_families/
├── FAM001/
│   ├── bgc_1.gbk
│   └── bgc_2.gbk
└── FAM002/
    └── ...
```

CDS features need `locus_tag`, `translation` and `gene_kind` qualifiers. Genes without `gene_kind` count as neither regulatory nor biosynthetic.

In the paper, families are gene cluster families (GCFs) from BiG-FAM, and each family directory holds the GenBank files of the GCF's member BGCs. `bin/group_gcf_families.py` builds this layout from a table of GCF ids and file names:

```
python bin/group_gcf_families.py gcf_table.tsv all_gbk_files/ bgc_families/
```

`gcf_table.tsv` is tab-separated with the columns `gcf_id` and `new_name` (GenBank file name). One `FAM<gcf_id>` directory is created per GCF.

## Running

```
git clone https://github.com/sysbiomelab/Motif-Miner.git
cd Motif-Miner

# locally with Singularity
nextflow run main.nf -profile standard --input_dir /path/to/bgc_families --outdir results

# locally with Docker
nextflow run main.nf -profile standard,docker --input_dir /path/to/bgc_families --outdir results

# SLURM cluster
nextflow run main.nf -profile hpc --input_dir /path/to/bgc_families --outdir results -resume
```

The `hpc` profile submits to the SLURM queue `cpu` with the resource allocations in `nextflow.config`; edit these for your cluster. On a shared cluster, set `NXF_CONDA_CACHEDIR` and `NXF_SINGULARITY_CACHEDIR` to a shared filesystem before launching. `nextflow run main.nf --help` lists all options.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | required | directory of `FAM*` subdirectories |
| `--outdir` | `./results` | final results |
| `--tempdir` | `./temp` | intermediate files from every step |
| `--min_upstream_len` | 100 | discard extracted regions shorter than this (bp) |
| `--max_upstream_len` | 350 | maximum upstream length to extract (bp) |
| `--inside_gene_len` | 20 | bases of the gene start included in each region |
| `--org_threshold` | 1.0 | a COG is assigned to a gene kind category when its percentage of that kind is above this |

### Fixed settings

These are hardcoded and were used as-is in the paper:

- Partial CDS features are included when converting GenBank to FASTA.
- ProteinOrtho runs with default parameters.
- CD-HIT-EST: `-c 0.90 -n 5 -d 0`.
- MEME runs only on COGs with more than 4 sequences after CD-HIT, with `-dna -revcomp -nmotifs 10`, in the six combinations of `-mod anr`/`-mod zoops` and `-maxw 10`/`20`/`30`.
- A COG is placed in a category when its percentage is strictly greater than `org_threshold`, and can be placed in more than one. COGs meeting no threshold go to `unknown`.
- PWM files use a uniform 0.25 background.

## Output

```
results/
├── families_excluded_no_meme_results.txt    only present if some families had no motifs
└── FAM001/
    ├── FAM001_complete_summary.tsv          every motif site with its gene annotation
    └── organised_results/
        ├── regulatory/
        │   └── COG_00005_upstream_350bp_inside_20bp_min_100bp_anr_maxw20_12seqs/
        │       ├── meme.xml, meme.txt, meme.html
        │       ├── sequences_info.csv
        │       ├── motifs_info.csv
        │       ├── combined_sequence_motif_binding.csv
        │       ├── <dir>_motif_1.csv
        │       └── <dir>_motif_1.pwm
        ├── biosynthetic/
        ├── biosynthetic_additional/
        └── unknown/
```

Each analysed COG has six result directories, named `<COG file>_<mode>_maxw<width>_<n>seqs`. Directories under `organised_results/` are copies, not symlinks.

Intermediates are kept under `<tempdir>/<FAM>/` in numbered folders: `1_cds_summary`, `2_faa_files`, `3_proteinortho_raw`, `3_proteinortho_cogs`, `4_upstream_filtered`, `5_cdhit`, `6_gene_info`, `7_meme_analysis`, `8_meme_analysis_parsed`.

## Citation

If you use this pipeline, please cite:

<!-- paper reference -->

## License

MIT. See [LICENSE](LICENSE).

## Contact

Idris Matine, idris.matine@kcl.ac.uk
