# MAD2026

Reproducible literature search accompanying the review article

> **Computational pathology approaches to tumour microenvironment characterisation in colorectal cancer**
> Vahid Anari, Eva Budinská, Vlad Popovici
> RECETOX, Faculty of Science, Masaryk University, Brno, Czech Republic
> *Molecular Aspects of Medicine* (2026), special issue *Data-Driven Diagnostics: Digital Pathology in the Age of AI*

---

## What this repository contains

A single script, `fetch_references.py`, that regenerates the candidate corpus
from which the studies discussed in the review were selected. It documents the
search so that the corpus can be inspected, criticised and rebuilt
independently, rather than being taken on trust.

The script reproduces everything up to, but not including, the manual
full-text selection, which is described in the article itself.

## Method

The literature relevant to this review is produced by several research
communities that describe overlapping work in almost disjoint vocabularies. No
single compact query spans them without an unmanageable loss of specificity, so
four targeted PubMed queries ("strands") are used and their union is taken.

| Strand | Scope | Why it is needed |
|---|---|---|
| **S1** | Computational pathology applied to the colorectal tumour microenvironment | The core literature of the review |
| **S2** | Quantitative digital assessment of established colorectal tissue biomarkers | This work is described either as automated quantification or as artificial intelligence, rarely both |
| **S3** | Pan-cancer spatial immune profiling | Methods validated partly on colorectal cohorts may never name the disease |
| **S4** | Inference of spatially resolved molecular information from histology | A small, technically distinct literature in bioinformatics venues, not anchored to a disease term |

The pipeline runs in three stages:

1. **Search** — the four strands are combined with `OR`; PubMed deduplicates the union.
2. **Harvest** — every record in the union is downloaded (PMID, year, journal, title, abstract, publication types, DOI).
3. **Screen** — deterministic rules are applied locally to titles and abstracts to remove records outside the scope of the review.

A record is removed at stage 3 if **at least one** rule rejects it; the removed
set is the union of the per-rule rejection sets, and every rule that rejected a
record is recorded in `excluded_records.csv`. Several rules fire only in the
absence of tissue-imaging vocabulary: studies included in the review do discuss
therapy, surgery and gene signatures, but all of them also examine tissue.

All queries and all screening rules are stated in full in the script, and the
full query strings are reproduced in the generated `search_summary.txt`.

## Installation

Requires Python 3.9 or later.

```bash
git clone https://github.com/bioinfo-recetox/MAD2026.git
cd MAD2026
pip install -r requirements.txt
```

## Usage

The NCBI usage policy requires a contact address. An API key is optional and
raises the request rate limit from 3/s to 10/s, which makes the download
noticeably faster.

```bash
export NCBI_EMAIL="your.address@institution.org"
export NCBI_API_KEY="..."        # optional

python3 fetch_references.py
```

Options:

| Option | Effect |
|---|---|
| `--outdir DIR` | Write results to `DIR` (default: `./search_YYYY-MM-DD`) |
| `--no-cache` | Re-download even if a harvest already exists |
| `--quiet` | Suppress progress bars (for logs and batch jobs) |
| `--selftest` | Run the offline tests of the screening logic and exit |

`--selftest` needs no network access and no API key, and verifies the screening
rules against a small set of constructed records in about a second.

The full run takes a few minutes, almost all of it in the download, which is
bounded by the NCBI rate limit. Progress is reported for both the download and
the screening stage.

## Outputs

| File | Contents |
|---|---|
| `candidate_records.csv` | **Main output** — the candidate list after screening |
| `candidate_records.txt` | The same list, formatted for reading |
| `excluded_records.csv` | Every removed record, with the rule(s) that removed it |
| `union_records.csv` | The complete unfiltered harvest (stage 2) |
| `search_summary.txt` | PRISMA-style counts, per-rule statistics, full query strings |
| `run_log.txt` | Timestamped log of the run |

## Reproducibility

PubMed is updated continuously, so the number of records retrieved depends on
the date of execution. The counts obtained by the authors are recorded in the
`REFERENCE_RUN` constant at the top of the script and printed in
`search_summary.txt` alongside the counts of the current run.

The queries, the screening rules and the procedure are fixed. A difference in
counts between runs therefore reflects growth of the database, not a change of
method.

A small number of studies discussed in the review are not retrievable by any
PubMed query — most obviously preprints that are not indexed there. These were
identified through citation tracking and expert knowledge, and are reported as
such in the article.

## Citation

If you use this code or the search strategy, please cite the article:

```bibtex
@article{anari2026computational,
  author  = {Anari, Vahid and Budinsk\'a, Eva and Popovici, Vlad},
  title   = {Computational pathology approaches to tumour microenvironment
             characterisation in colorectal cancer},
  journal = {Molecular Aspects of Medicine},
  year    = {2026},
  note    = {Special issue: Data-Driven Diagnostics: Digital Pathology in the Age of AI}
}
```

## License

MIT — see [LICENSE](LICENSE).

## Contact

Vlad Popovici — RECETOX, Faculty of Science, Masaryk University, Brno, Czech Republic
