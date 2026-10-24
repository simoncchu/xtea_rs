# xtea-rs: xTea in Rust (Illumina short reads)

`xtea-rs` re-implements the **master branch** of [xTea](https://github.com/parklab/xTea)
(comprehensive transposable element analyzer) in Rust, for **paired-end Illumina short reads**
(WGS or WES). The binaries are drop-in replacements for the Python programs. They take the same
options, write the same pipeline layout and produce the same output files. With the same random
seed, every result file compared so far is byte-identical to the Python version
(see [Validation](#validation)).

Python is not needed at run time. The DF21 genotyping model (a joblib pickle) is read by a
built-in pickle reader, and the deep-forest inference is re-implemented.

## Scope

**Ported** (everything that the master branch runs on Illumina data):

* the pipeline generators for hg38, hg19 and CHM13, plus the single-sample cloud wrapper;
* every step of `x_TEA_main.py` that the generated pipelines run:
  * clipped reads (`-C`) and discordant pairs (`-D`);
  * consensus filtering (`-N`);
  * 3' transductions (`--transduction`) and sibling/orphan transductions (`--sibling`);
  * post-filtering (`--postF`, `--postFmosaic`);
  * gene annotation (`--gene`);
  * DF21 genotyping (`--gntp_classify`, `--gntp_feature`);
  * VCF output (`--gVCF`);
* the other modes of `x_TEA_main.py`:
  * case-control/somatic calling (`--case_control`, `--somatic_hc`);
  * mosaic calling (`--mosaic`) and joint mosaic filtering (`-J/--joint`);
  * IGV scripts (`--igv`);
  * library preprocessing (`-P`, `--flank`);
  * clipped-read collection (`-Q`);
* `x_TEA_main_germline_cloud.py`;
* the stand-alone tools `x_vcf_merger.py`, `x_cmp.py` and `x_population_AF.py`.

**Not ported:**

* **10X linked reads.** The Rust version supports Illumina short reads only:
  * the 10X modes of `xtea_main` exit with an error. These are `-B/--barcode`,
    `-E/--collect`, `-I/--mutation`, `-A/--assembly`, `-M/--map`, `--flk_map` and
    `-F/--filter_asm`;
  * a BAM-list entry tagged `10X` is rejected by the clip step;
  * the generators reject an existing `-x` 10X BAM list and the 10X steps of `-f`
    (4, 8, 32, 64 and 128).

  Use the Python xTea for 10X data.
* **Long reads.** `bin/xtea_long` generates commands for the separate xTea_long project.
* **Training the genotype classifier** (`--train_gntp`). The original training code depends on
  site lists from the authors' cluster and on scikit-learn. Prediction with the shipped DF21
  model is fully supported.

A few options exist in `x_TEA_main.py` but have no code behind them, such as `--sv`,
`--spectrum`, `--analysis` and `-U`. They are accepted and do nothing, as in the original.

## Building

```bash
cd xtea-rs
cargo build --release                     # binaries in target/release/
cargo build --release --features remote   # optional: http/https/s3 BAM/CRAM via htslib + libcurl
cargo test --release                      # unit and parity tests
```

The build needs:

* a recent stable Rust toolchain (tested with 1.97);
* a C compiler, `make` and CMake. htslib, zlib-ng, bzip2 and xz are compiled from source and
  linked statically, so no system htslib is needed.

It has been built and tested on Linux x86_64. macOS should work with the Xcode command-line tools
and CMake (`brew install cmake`), but has not been tested.

At run time you need:

* `bwa` and `samtools` on `PATH`, as for the Python version;
* the repeat library (`rep_lib_annotation`), the reference genome and the gene annotation;
* the `xtea/` folder of this repository, passed with `--xtea` as before. It only provides the
  genotyping model `genotyping/DF21_model_1_2`.

## Programs

| Rust binary | Python original | Purpose |
|---|---|---|
| `xtea` | `bin/xtea` = `xtea/gnrt_pipeline_local_v38.py` | pipeline generator, hg38 |
| `xtea_hg19` | `bin/xtea_hg19` = `xtea/gnrt_pipeline_local.py` | pipeline generator, hg19 |
| `xtea_chm13` | `xtea/gnrt_pipeline_local_chm13.py` | pipeline generator, CHM13 |
| `xtea_germline_cloud` | `xtea/gnrt_pipeline_germline_cloud_v38.py` | single-sample germline wrapper of the cloud workflow |
| `xtea_main` | `xtea/x_TEA_main.py` | the calling steps run by the generated pipelines |
| `xtea_main_germline_cloud` | `xtea/x_TEA_main_germline_cloud.py` | the cloud variant of `xtea_main` |
| `xtea_vcf_merger` | `xtea/x_vcf_merger.py` | merge per-sample VCFs |
| `xtea_cmp` | `xtea/x_cmp.py` | compare call sets |
| `xtea_population_af` | `xtea/x_population_AF.py` | population allele frequencies from a merged VCF |

Keep the binaries together in one folder (for example `target/release/`, or copy them all to one
`bin` folder). The generators write the path of the `xtea_main` that sits next to them into the
generated scripts. If none is found there, they fall back to `xtea_main` on `PATH`.

## Usage

Use the Rust binaries the same way as the Python ones. For example, a germline run on hg38 for
L1, Alu and SVA:

```bash
xtea -i sample_id.txt -b illumina_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh \
     -l /path/rep_lib_annotation/ -r /path/reference/genome.fa -g /path/gencode.gff3 \
     --xtea /path/xTea/xtea/ -f 5907 -y 7 --slurm -t 0-12:00 -q short -n 8 -m 25
sh submit_jobs.sh    # or submit the generated run_xTEA_pipeline.sh scripts yourself
```

Case-control (somatic) mode on tumour/normal pairs:

```bash
xtea --case_control --tumor -i sample_id.txt -b case_ctrl_bam_list.txt -p ./path_work_folder/ \
     -o submit_jobs.sh -l /path/rep_lib_annotation/ -r /path/reference/genome.fa \
     -g /path/gencode.gff3 --xtea /path/xTea/xtea/ -f 5907 -y 7 --slurm -n 8 -m 25
```

Input files, flags and outputs are as described in the main [xTea README](../README.md), with
these differences in the generated pipelines:

* Each command calls `${XTEA_MAIN}` instead of `python ${XTEA_PATH}"x_TEA_main.py"`. Every
  `<rep>.config` gains a line `XTEA_MAIN <path>` right after `XTEA_PATH`.
* `-x` (the 10X BAM list) is optional and only `null` is accepted. The configurations no longer
  define `BAM1` and `BARCODE_BAM`.
* The `-I` (10X internal-mutation) step is not generated. Without this change, flag 256 would add
  it to the HERV, MSTA and Mitochondrion scripts.
* `-f`: 1 clip, 2 discordant, 16 filter, 256 transduction, 512 post-filter, 1024 annotation and
  genotyping (+2048 skips it), 4096 sibling transduction. `5907` runs all of them. As in the
  Python version, `-f` has no default in `xtea` and `xtea_chm13`. In `xtea_hg19` it defaults
  to 19.
* `-y`: 1 L1, 2 Alu, 4 SVA, 8 HERV, 16 Mitochondrion, 32 MSTA, 64 Pseudogene (hg38/CHM13 only).
* `-V` prints `xTea v0.1 for Illumina short reads on <genome>`.

`xtea_main -h` lists all steps and options. Help texts of the 10X-only options are marked
`[10X; not supported]` or `[10X; unused]`. Usage errors are reported like Python's `optparse`
(`xtea_main: error: no such option: --x`, exit status 2), including unique-prefix matching of long
options.

### Random seed and reproducibility

Before the clip, discordant and filter steps, xTea estimates coverage, read length and insert size
from randomly sampled sites. Unless `--user` is given, these estimates set the cutoffs. Python
samples the sites with an unseeded `random`, so runs differ slightly from one another. The Rust
binaries do the same by default. With `xtea_main --seed N`, or `XTEA_SEED=N` in the environment,
they draw exactly the sequence CPython draws after `random.seed(N)`. A seeded run is therefore
reproducible, and it can be compared one-to-one with a Python run seeded the same way. That is how
the validation below was done.

## Python behaviour kept on purpose

The port reproduces the Python master's behaviour, including its quirks and failure modes, so the
two versions can be compared file by file. The ones users are likely to meet:

* **BamSnap steps fail (hg38/CHM13).** The annotation stage ends with two BamSnap screenshot
  commands, which stop with `no such option: --bamsnap`. `x_TEA_main.py` has no such option. The
  calls and VCFs are produced before these commands and are not affected.
* **Options the master `x_TEA_main.py` does not have.** The following generator settings emit
  commands that fail with `no such option`:
  * `--hard` adds `--hard` to the post-filter commands (only `x_TEA_main_germline_cloud.py`
    accepts it);
  * flag 8192 emits `--sr_asm`;
  * de novo mode emits `--denovo`.
* **Mosaic post-filter needs a blacklist.** `--postF --postFmosaic` filters the calls against a
  blacklist. With the default `null` it stops because the file `null` cannot be opened (in
  Python, a `FileNotFoundError`). If a chromosome has calls but is missing from the blacklist,
  it stops with the original `TypeError: write() argument must be str, not int`.
* **Joint mosaic filtering (`--joint`).** Two of the original membership tests compare a tuple
  instead of a boolean, which changes which singletons are reported. The port keeps this; the
  details are in `src/mosaic.rs`.
* **Python-style messages.** Errors are reported the way Python raises them, for example
  `TypeError` when `-f` is missing, or `IndexError` and `ValueError` on malformed input lines.
* **Python semantics throughout.** Number formatting, `int()`/`float()` parsing, `str.split()`
  whitespace, universal newlines, set/dict iteration order and true division all follow
  CPython. This is what makes the outputs byte-identical.

## Validation

All comparisons were made against the Python master (Python 3.8). Unseeded randomness was seeded
on both sides with the same value.

The test data is a simulated paired-end library: about 2.0 M aligned reads on a 10 Mb reference
(chr1, chr2, chrM and a decoy contig). It contains 39 insertions: 14 L1, 2 orphan L1
transductions, 16 Alu and 7 SVA.

| What | Result |
|---|---|
| Germline pipeline, L1/Alu/SVA, all 10 steps per type (clip → gVCF), scripts made by the Rust generator | 202 of 202 result files identical (only the config and script files differ, by design) |
| Case-control mode (the test sample against a simulated control) | identical |
| Mosaic post-filter (`--postFmosaic`) and `--joint` | identical |
| DF21 genotype classifier, 20,000 feature vectors | bit-identical probabilities and genotypes |
| Pipeline generators, 3,450 option combinations | byte-identical scripts and configs (after the Illumina-only change, they differ only in the removed 10X lines) |
| `x_vcf_merger` / `x_cmp` / `x_population_AF`, 133 cases | identical (fixtures in `tests/data/tools`) |

On this small data set the three pipelines took 33 s, 19 s and 20 s with the Rust binaries,
against 104 s, 108 s and 82 s with Python (`-n 4`, same machine). Much of the remaining time is
spent in `bwa` and `samtools`. Whole-genome runs have not been benchmarked yet.

## Code layout

| Path | Contents |
|---|---|
| `src/bin/` | one file per binary |
| `src/cli.rs`, `src/steps.rs` | `xtea_main` options and the per-mode drivers (the `if/elif` chain of `x_TEA_main.py`) |
| `src/pipeline/` | the four pipeline generators |
| `src/tools/` | the stand-alone tools |
| `src/*.rs` | one module per Python module, e.g. `tei_locator.rs` ← `x_TEI_locator.py`, `post_filter.rs` ← `x_post_filter.py`, `genotype_classify.rs` ← `x_genotype_classify.py` + DF21 inference |
| `src/bam.rs`, `src/fasta.rs` | thin pysam-like layer over rust-htslib |
| `src/py/` | CPython/numpy semantics: set hashing and order, `random` (MT19937), float formatting, `optparse`, pickle, text handling |
| `tests/` | parity tests (CPython set order, intervaltree, tools fixtures) |

## License

Same terms as xTea: see [`LICENSE`](../LICENSE) (Academic and Non-Commercial Research Use
Software License).
