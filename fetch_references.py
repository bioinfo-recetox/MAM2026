#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_references.py
===============================================================================
Reproducible literature search for:

    "Computational pathology approaches to tumour microenvironment
     characterisation in colorectal cancer"
    V. Anari, E. Budinska, V. Popovici
    RECETOX, Faculty of Science, Masaryk University, Brno, Czech Republic

WHAT THIS SCRIPT DOES
---------------------
It regenerates the candidate corpus from which the studies reviewed in the
article were selected, in three stages:

  Stage 1  SEARCH    Four complementary PubMed queries ("strands"), each
                     targeting a research community that publishes on this
                     topic in its own vocabulary, are combined with OR.
                     PubMed deduplicates the union.
  Stage 2  HARVEST   Every record in the union is downloaded (PMID, year,
                     journal, title, abstract, publication types, DOI).
  Stage 3  SCREEN    Deterministic rules are applied locally to titles and
                     abstracts, removing records outside the scope of the
                     review.

The output of Stage 3 is the candidate list that the three authors then
screened manually, in full text, to select the studies discussed in the
review. This script reproduces everything up to, but not including, that
manual selection, which is described in the article itself.

USAGE
-----
    pip install requests            # tqdm is optional, see below
    export NCBI_EMAIL="your.address@institution.org"   # required by NCBI
    export NCBI_API_KEY="..."                          # optional, faster
    python3 fetch_references.py

  Options
    --outdir DIR     write results to DIR (default: ./search_YYYY-MM-DD)
    --no-cache       re-download even if a harvest already exists
    --quiet          suppress progress bars (for logs and batch jobs)
    --selftest       run the offline tests of the screening logic and exit

  The run takes a few minutes: the harvest is limited by NCBI's rate limit,
  and the screening rules are applied to every harvested abstract. Progress
  is reported for both stages. If the 'tqdm' package is installed it is used
  for the progress bars; otherwise an equivalent built-in bar is shown.

OUTPUTS (in --outdir)
---------------------
    candidate_records.csv   final candidate list after screening  <- MAIN OUTPUT
    candidate_records.txt   the same list, formatted for reading
    excluded_records.csv    every removed record, with the rule(s) that removed it
    union_records.csv       the complete unfiltered harvest (Stage 2)
    search_summary.txt      PRISMA-style counts, per-rule statistics, full queries
    run_log.txt             timestamped log of the run

REPRODUCIBILITY NOTE
--------------------
PubMed is updated continuously, so the number of records retrieved depends on
the date of execution. The counts obtained by the authors are recorded in
REFERENCE_RUN below. The queries, the rules and the procedure are fixed, so a
difference in counts reflects database growth rather than a change of method.

This script is released with the article so that the corpus can be inspected,
criticised and rebuilt independently.
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import re
import sys
import time
import xml.etree.ElementTree as ET

try:
    import requests
except ImportError:  # pragma: no cover
    sys.exit("This script requires the 'requests' package:  pip install requests")

try:
    from tqdm import tqdm
    _HAVE_TQDM = True
except ImportError:
    _HAVE_TQDM = False


# =============================================================================
# REFERENCE RUN -- counts obtained by the authors
# =============================================================================

REFERENCE_RUN = {
    "date": "2026-09-15",
    "union_records": 9228,        
    "candidates_after_screening": 785,
    "studies_selected_manually": 49,
}


###############################################################################
#  D O   N O T   C H A N G E   T H E   R U L E S   B E L O W   U N L E S S    
#  Y O U   K N O W   W H A T   Y O U   D O!!!
###############################################################################

# =============================================================================
# STAGE 1 -- SEARCH STRATEGY
# =============================================================================
# The relevant literature is produced by several research communities that
# describe overlapping work in almost disjoint vocabularies. No single compact
# query spans them without an unmanageable loss of specificity, so four
# targeted strands are used and their union is taken.
# =============================================================================

DATE_RANGE = '("2016/01/01"[dp] : "2026/12/31"[dp])'
LANGUAGE = "english[la]"

_CRC = (
    '(colorect*[tiab] OR colon*[tiab] OR rect*[tiab] OR CRC[tiab] '
    'OR "Colorectal Neoplasms"[MeSH])'
)

# --- Strand 1 ---------------------------------------------------------------
# Computational pathology applied to the colorectal tumour microenvironment:
# the core literature of the review. The tissue block includes cell-marker and
# immunohistochemistry terms, because studies from the diagnostic-pathology
# tradition describe their target as CD3/CD8 T cells rather than as "immune
# infiltrate" (e.g. the QUASAR analysis, PMID 39083705).
STRAND_1 = (
    _CRC + " AND "
    '(microenviron*[tiab] OR immune[tiab] OR lymphocyt*[tiab] OR lymphoid[tiab] '
    'OR "lymph node"[tiab] OR "lymph nodes"[tiab] OR strom*[tiab] '
    'OR fibroblast*[tiab] OR collagen[tiab] OR budding[tiab] OR eosinophil*[tiab] '
    'OR CD3[tiab] OR CD8[tiab] OR FOXP3[tiab] OR immunohistochem*[tiab] '
    'OR "T-cell"[tiab] OR "T cells"[tiab] OR "tumor infiltrating"[tiab] '
    'OR "tumour infiltrating"[tiab])'
    " AND "
    '("deep learning"[tiab] OR "machine learning"[tiab] '
    'OR "artificial intelligence"[tiab] OR AI[tiab] '
    'OR "computational pathology"[tiab] OR "digital pathology"[tiab] '
    'OR "whole slide"[tiab] OR neural[tiab])'
    f" AND {DATE_RANGE} AND {LANGUAGE}"
)

# --- Strand 2 ---------------------------------------------------------------
# Quantitative digital assessment of established colorectal tissue biomarkers.
# This work is described either in the vocabulary of automated quantification
# or in that of artificial intelligence, rarely both, so the method clause
# accepts either. Records relying solely on rule-based thresholding, template
# counting or manual stereology are removed at full-text screening, where the
# eligibility criterion for AI-based modelling is applied.
STRAND_2 = (
    _CRC + " AND "
    '(immunoscore[tiab] OR "tumor-stroma ratio"[tiab] OR "tumour-stroma ratio"[tiab] '
    'OR "tumor infiltrating"[tiab] OR "tumour infiltrating"[tiab] OR CD3[tiab] '
    'OR CD8[tiab] OR FOXP3[tiab] OR "PD-1"[tiab] OR immunohistochem*[tiab] '
    'OR "tumor budding"[tiab] OR "tumour budding"[tiab] OR collagen[tiab] '
    'OR "lymph node"[tiab] OR "lymph nodes"[tiab])'
    " AND "
    '(("image analysis"[tiab] OR "digital image"[tiab] OR "digital pathology"[tiab] '
    'OR automat*[tiab] OR quantitat*[tiab] OR quantif*[tiab] OR algorithm*[tiab] '
    'OR stereolog*[tiab] OR "computer-assisted"[tiab]) '
    'OR ("deep learning"[tiab] OR "machine learning"[tiab] '
    'OR "artificial intelligence"[tiab] OR AI[tiab] OR neural[tiab] '
    'OR "computational pathology"[tiab] OR convolutional[tiab] '
    'OR "self-supervised"[tiab] OR "weakly supervised"[tiab]))'
    f" AND {DATE_RANGE} AND {LANGUAGE}"
)

# --- Strand 3 ---------------------------------------------------------------
# Pan-cancer spatial immune profiling. Methods developed across tumour types
# and validated partly on colorectal cohorts may never name the disease in
# title or abstract; their relevance is established at screening.
STRAND_3 = (
    '("pan-cancer"[tiab] OR pancancer[tiab] OR "solid tumors"[tiab] '
    'OR "solid tumours"[tiab] OR "cancer types"[tiab] '
    'OR "multiple cancers"[tiab])'
    " AND "
    '(CD8[tiab] OR "PD-1"[tiab] OR "tumor infiltrating"[tiab] '
    'OR "tumour infiltrating"[tiab] OR lymphocyt*[tiab] OR immune[tiab])'
    " AND "
    '("whole slide"[tiab] OR "deep learning"[tiab] OR "machine learning"[tiab] '
    'OR "artificial intelligence"[tiab] OR "image analysis"[tiab] '
    'OR "digital pathology"[tiab] OR quantitat*[tiab] OR immunoprofil*[tiab] '
    'OR "spatial heterogeneity"[tiab] OR "spatial patterns"[tiab])'
    f" AND {DATE_RANGE} AND {LANGUAGE}"
)

# --- Strand 4 ---------------------------------------------------------------
# Inference of spatially resolved molecular information from histology. A small
# and technically distinct literature, published largely in bioinformatics
# venues, deliberately not anchored to a disease term.
STRAND_4 = (
    '("spatial transcriptomics"[tiab] OR "spatially resolved transcriptomics"[tiab])'
    " AND "
    '(histolog*[tiab] OR histopatholog*[tiab] OR "H&E"[tiab] OR "whole slide"[tiab])'
    " AND "
    '("deep learning"[tiab] OR "machine learning"[tiab] OR "cross-modal"[tiab] '
    'OR "contrastive learning"[tiab] OR neural[tiab] OR predict*[tiab])'
    f" AND {DATE_RANGE} AND {LANGUAGE}"
)

STRANDS = {
    "S1_crc_tme_computational_pathology": STRAND_1,
    "S2_crc_biomarker_quantification": STRAND_2,
    "S3_pancancer_spatial_immune": STRAND_3,
    "S4_spatial_omics_from_histology": STRAND_4,
}

# The outer parentheses are required: PubMed binds AND more tightly than OR,
# so an unwrapped disjunction changes meaning when combined with other terms.
UNION_QUERY = "((" + ") OR (".join(STRANDS.values()) + "))"


# =============================================================================
# STAGE 3 -- SCREENING RULES
# =============================================================================
# Applied to the title and abstract of every harvested record, case-insensitive.
#
#   EXCLUDE_RULES  remove a record when the pattern IS found.
#   REQUIRE_RULES  remove a record when the pattern is NOT found, within the
#                  stated scope ("abstract" or "title_abstract"). A record with
#                  no abstract is therefore removed by abstract-scoped rules,
#                  since it cannot be screened on its content.
#
# A record is removed if at least one rule rejects it; the set of removed
# records is the union of the per-rule rejection sets, and every rule that
# rejected a record is recorded in excluded_records.csv.
#
# Several rules use a negative lookahead for tissue-imaging vocabulary, e.g.
#   (?is)^(?!.*\b<imaging terms>\b).*\b<topic terms>\b
# These fire only when the record shows no evidence of examining tissue. The
# guard is essential: studies included in the review do discuss therapy,
# surgery and gene signatures, but all of them also image tissue.
# =============================================================================

EXCLUDE_RULES = [
    (
        "other_imaging_modality",
        r"\b(radiomic\w*|computed tomograph\w*|magnetic resonance imag\w*|MRI"
        r"|CT scan\w*|ultrasound|ultrasonograph\w*|positron emission|x-ray"
        r"|mass[ -]spectrometry|western blot|immunoblot)\b",
    ),
    (
        "endoscopy_or_screening",
        r"\b(colonoscop\w+|endoscop\w+|polyp detection|capsule endoscopy"
        r"|f(?:e|ae)cal|stool|screening program\w*|precancer\w*|laparoscop\w+)\b",
    ),
    (
        "liquid_biopsy_or_circulating",
        r"\b(liquid biops\w+|ctDNA|cell-free DNA|circulating tumou?r"
        r"|peripheral blood|serum|plasma(?!\s*cell)|lncRNA\w*|miRNA\w*)\b",
    ),
    (
        "quality_of_life_or_survey",
        r"\b(quality of life|questionnaire\w*)\b",
    ),
    (
        "preclinical_models",
        r"\b(xenograft\w*|organoid\w*|zebrafish|murine|mouse|mice|rat|rats"
        r"|spheroid\w*|canine|porcine|in vivo|ex vivo|in vitro|fung\w+)\b",
    ),
    (
        "other_organ_or_disease",
        r"\b(brain|glioblastoma|gleason|ampullary|lung adenocarcinoma|renal"
        r"|pituitary|gastric cancer|pancreatic cancer|breast cancer"
        r"|Lynch syndrome|nerve)\b",
    ),
    (
        "other_molecular_topic",
        r"\b(senescence|ageing|haematolog\w+|single[ -]cell|developmental"
        r"|pathomics|polypeptide|peptide\w*)\b",
    ),
    # Another primary tumour is the subject and colorectal cancer is never
    # mentioned. The leading negative lookahead protects pan-cancer studies
    # that include a colorectal cohort, which are in scope.
    (
        "other_primary_tumour",
        r"(?is)^(?!.*\b(colorect\w*|colon|colonic|rectal|rectum|crc)\b)"
        r".*\b(breast cancer|non.small.cell lung|nsclc|lung adenocarcinom\w*"
        r"|prostat\w+ cancer|gastric cancer|hepatocellular"
        r"|pancreatic (ductal )?adenocarcinom\w*|ovarian cancer|melanoma"
        r"|head and neck (squamous|cancer)|bladder cancer|urothelial"
        r"|renal cell carcinom\w*|cervical cancer|e?sophageal cancer|glioma"
        r"|glioblastoma|lymphoma|leukaemia|leukemia)\b",
    ),
    # Review, survey or protocol language that the publication-type filter
    # misses, e.g. preprints and records typed only as "Journal Article".
    (
        "review_or_survey_language",
        r"\b(this review|we review|systematic review|meta.analysis"
        r"|narrative review|a survey|study protocol|protocol for a)\b",
    ),
    # A wet-lab molecular assay is the method of the study.
    (
        "wet_lab_assay",
        r"\b(reverse transcription|rt-?pcr|qpcr|real-?time pcr"
        r"|polymerase chain reaction|western blot\w*|immunoblot\w*|elisa"
        r"|flow cytometr\w+|cell line\w*|transfect\w+|knockdown|knock-?out"
        r"|sirna|shrna|luciferase|co-?culture|cck-?8|wound healing assay"
        r"|transwell|colony formation)\b",
    ),
]

REQUIRE_RULES = [
    # The study must involve observation of tissue. A bare "imag*" pattern is
    # too narrow: papers from the diagnostic-pathology tradition describe their
    # material as immunohistochemistry, staining or tissue sections and may
    # never use an imaging word.
    (
        "tissue_observation_required",
        r"\b(imag\w*|histolog\w*|histopatholog\w*|immunohistochem\w*"
        r"|immunofluorescen\w*|stain\w*|microscop\w*|whole[- ]slide|H&E"
        r"|h and e|digital patholog\w*|tissue section\w*|biops\w*"
        r"|morpholog\w*|slide\w*|specimen\w*|immunoscore|immune context\w*)\b",
        "title_abstract",
    ),
    # Some computational or quantitative method must be present. Deliberately
    # broad, because it must admit the Immunoscore and stereology literature,
    # which speaks of densities and scores rather than of learning.
    (
        "method_vocab_required",
        r"\b(deep learning|machine learning|artificial intelligence|ai"
        r"|neural network\w*|convolutional|cnn|transformer\w*|computational|comput\w+"
        r"|digital patholog\w*|whole[- ]slide|wsi|image analysis|automat\w+|algorithm\w*"
        r"|segmentat\w+|classif\w+|self[- ]supervised|weakly[- ]supervised"
        r"|foundation model\w*|graph neural|random forest|support vector|regression"
        r"|predictive model\w*|quantif\w+|quantitat\w+|densit\w+|scoring|score\w*"
        r"|stereolog\w+|data[- ]driven|software|model\w*|learning)\b",
        "title_abstract",
    ),
    # The study must concern the tumour microenvironment, not merely tissue.
    # Without this, molecular-marker, staging and surgical-outcome studies
    # survive on the strength of a single mention of immunohistochemistry.
    (
        "tme_vocab_required",
        r"\b(microenviron\w*|immune|immunolog\w*|lymphocyt\w*|til|tils"
        r"|stroma\w*|fibroblast\w*|caf|cafs|collagen|extracellular matrix|ecm"
        r"|budding|cd3|cd8|cd4|foxp3|pd-?l?1|macrophag\w*|neutrophil\w*|eosinophil\w*"
        r"|immunoscore|infiltrat\w+|tumou?r[- ]stroma|invasive margin|invasive front"
        r"|desmoplas\w+|niche\w*|neighbou?rhood\w*|cell type\w*|spatial"
        r"|protein expression|tissue architect\w*|contexture)\b",
        "title_abstract",
    ),
    # Colorectal cancer must be named, except for disease-agnostic methods that
    # infer spatial molecular information from histology; those papers never
    # name the disease in title or abstract, and their relevance is established
    # at full text.
    (
        "crc_or_st_method_required",
        r"\b(colorect\w*|colon|colonic|rectal|rectum|crc|coad|read|pan[- ]cancer"
        r"|pancancer|across \d+ cancers|multiple cancer\w*|solid tumou?rs?"
        r"|cancer types|tcga)\b"
        # The "^" anchor on this branch is load-bearing, not cosmetic. Without
        # it the engine retries the lookahead pair at every character position,
        # rescanning the whole abstract each time; that one omission accounted
        # for roughly 85% of the total screening time. [\s\S] is used rather
        # than "." so the branch does not depend on the DOTALL flag.
        r"|^(?=[\s\S]*\bspatial transcriptomic\w*\b)"
        r"(?=[\s\S]*\b(histolog\w*|histopatholog\w*|whole[- ]slide|wsi)\b)",
        "title_abstract",
    ),
]

EXCLUDED_PUBLICATION_TYPES = {
    "Review", "Systematic Review", "Meta-Analysis", "Editorial", "Comment",
    "Letter", "Published Erratum", "Retraction of Publication", "Case Reports",
    "Clinical Trial Protocol", "Guideline", "News",
}


# =============================================================================
# Implementation
# =============================================================================

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_RETRIES = 4
BATCH_SIZE = 200

_log_path: str | None = None
_quiet = False


def log(message: str) -> None:
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    line = f"[{stamp}] {message}"
    if not _quiet:
        print(line, flush=True)
    if _log_path:
        with open(_log_path, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


class Progress:
    """Progress bar: tqdm when available, otherwise a built-in equivalent.

    Written to stderr so that redirecting stdout to a file keeps the output
    clean. Silent when --quiet is given or when stderr is not a terminal.
    """

    def __init__(self, total: int, description: str):
        self.total = max(1, total)
        self.n = 0
        self.description = description
        self.start = time.time()
        self.enabled = not _quiet and sys.stderr.isatty()
        self.bar = None
        if self.enabled and _HAVE_TQDM:
            self.bar = tqdm(total=total, desc=description, unit="rec",
                            file=sys.stderr, leave=False)

    def update(self, step: int = 1) -> None:
        self.n += step
        if not self.enabled:
            return
        if self.bar is not None:
            self.bar.update(step)
            return
        if self.n % 25 and self.n != self.total:
            return
        frac = self.n / self.total
        width = 32
        filled = int(width * frac)
        elapsed = time.time() - self.start
        eta = (elapsed / frac - elapsed) if frac > 0 else 0.0
        sys.stderr.write(
            f"\r  {self.description:24s} |{'#' * filled}{'.' * (width - filled)}| "
            f"{self.n:>6d}/{self.total} ({frac:4.0%})  ETA {eta:4.0f}s"
        )
        sys.stderr.flush()

    def close(self) -> None:
        if not self.enabled:
            return
        if self.bar is not None:
            self.bar.close()
            return
        elapsed = time.time() - self.start
        sys.stderr.write(
            f"\r  {self.description:24s} |{'#' * 32}| "
            f"{self.n:>6d}/{self.total} (100%)  done in {elapsed:.0f}s\n"
        )
        sys.stderr.flush()


def _params(extra: dict) -> dict:
    p = {
        "tool": "crc_tme_review_fetch_references",
        "email": os.environ.get("NCBI_EMAIL", ""),
    }
    api_key = os.environ.get("NCBI_API_KEY", "")
    if api_key:
        p["api_key"] = api_key
    p.update(extra)
    return p


def _delay() -> float:
    return 0.12 if os.environ.get("NCBI_API_KEY") else 0.40


def _request(endpoint: str, params: dict, timeout: int = 180):
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(f"{EUTILS}/{endpoint}",
                                    params=_params(params), timeout=timeout)
            if response.status_code == 200:
                time.sleep(_delay())
                return response
            last_error = f"HTTP {response.status_code}"
        except Exception as exc:
            last_error = str(exc)
        time.sleep(_delay() * (2 ** (attempt + 1)))
    raise RuntimeError(f"NCBI request failed after {MAX_RETRIES} attempts ({last_error})")


def count_hits(query: str) -> int:
    response = _request("esearch.fcgi",
                        {"db": "pubmed", "term": query, "retmode": "json", "retmax": 0})
    return int(response.json()["esearchresult"]["count"])


def _text_of(element) -> str:
    if element is None:
        return ""
    return " ".join(ET.tostring(element, method="text", encoding="unicode").split())


def _parse_records(xml_bytes: bytes) -> list[dict]:
    """Parse journal articles and book chapters.

    PubMed returns book chapters as <PubmedBookArticle>, which esearch counts
    but which use a different element structure; ignoring them loses records.
    XPaths are anchored, because unanchored ".//PMID" and ".//ArticleId" also
    match the identifiers of cited references in PMC-deposited records.
    """
    root = ET.fromstring(xml_bytes)
    out = []

    for art in root.findall(".//PubmedArticle"):
        doi = ""
        for ident in art.findall("./PubmedData/ArticleIdList/ArticleId"):
            if ident.get("IdType") == "doi":
                doi = (ident.text or "").strip()
                break
        if not doi:
            for eloc in art.findall("./MedlineCitation/Article/ELocationID"):
                if eloc.get("EIdType") == "doi":
                    doi = (eloc.text or "").strip()
                    break
        out.append({
            "pmid": (art.findtext("./MedlineCitation/PMID") or "").strip(),
            "year": (art.findtext(".//Journal//PubDate/Year")
                     or art.findtext(".//Journal//PubDate/MedlineDate") or "")[:4],
            "journal": art.findtext(".//Journal/ISOAbbreviation")
                       or art.findtext(".//Journal/Title") or "",
            "title": _text_of(art.find(".//ArticleTitle")),
            "abstract": " ".join(_text_of(n) for n in
                                 art.findall(".//Abstract/AbstractText")).strip(),
            "publication_types": ";".join((n.text or "") for n in
                                          art.findall(".//PublicationType")),
            "doi": doi,
            "record_type": "journal",
        })

    for book in root.findall(".//PubmedBookArticle"):
        doi = ""
        for ident in book.findall("./PubmedData/ArticleIdList/ArticleId"):
            if ident.get("IdType") == "doi":
                doi = (ident.text or "").strip()
                break
        out.append({
            "pmid": (book.findtext("./BookDocument/PMID") or "").strip(),
            "year": (book.findtext(".//Book/PubDate/Year") or "")[:4],
            "journal": _text_of(book.find(".//Book/BookTitle")),
            "title": (_text_of(book.find(".//BookDocument/ArticleTitle"))
                      or _text_of(book.find(".//Book/BookTitle"))),
            "abstract": " ".join(_text_of(n) for n in
                                 book.findall(".//Abstract/AbstractText")).strip(),
            "publication_types": ";".join((n.text or "") for n in
                                          book.findall(".//PublicationType")),
            "doi": doi,
            "record_type": "book",
        })

    return out


def harvest_union() -> list[dict]:
    response = _request(
        "esearch.fcgi",
        {"db": "pubmed", "term": UNION_QUERY, "retmode": "json",
         "usehistory": "y", "retmax": 0},
    )
    result = response.json()["esearchresult"]
    total = int(result["count"])
    webenv, query_key = result["webenv"], result["querykey"]
    log(f"Union query returns {total} records; downloading")

    records: list[dict] = []
    progress = Progress(total, "downloading records")
    for start in range(0, total, BATCH_SIZE):
        response = _request(
            "efetch.fcgi",
            {"db": "pubmed", "WebEnv": webenv, "query_key": query_key,
             "retstart": start, "retmax": BATCH_SIZE, "retmode": "xml"},
        )
        parsed = _parse_records(response.content)
        records.extend(parsed)
        progress.update(min(BATCH_SIZE, total - start))
    progress.close()

    seen, unique = set(), []
    for record in records:
        if record["pmid"] and record["pmid"] not in seen:
            seen.add(record["pmid"])
            unique.append(record)

    if len(unique) != total:
        log(f"NOTE: {total} records reported, {len(unique)} parsed. PubMed "
            f"occasionally reports records that efetch cannot return "
            f"(withdrawn or in-process entries).")
    return unique


def _scope_text(record: dict, scope: str) -> str:
    if scope == "abstract":
        return record.get("abstract", "")
    return record.get("title", "") + " " + record.get("abstract", "")


def apply_rules(records: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    """Apply the screening rules.

    A record is removed if at least one rule rejects it, so the removed set is
    the union of the per-rule rejection sets. Returns (kept, excluded, stats).
    """
    exclude_compiled = [(name, re.compile(pattern, re.I))
                        for name, pattern in EXCLUDE_RULES]
    require_compiled = [(name, re.compile(pattern, re.I), scope)
                        for name, pattern, scope in REQUIRE_RULES]

    rejected_by_rule: dict[str, set] = {}
    kept, excluded = [], []

    progress = Progress(len(records), "screening records")
    for record in records:
        reasons = []
        for name, regex in exclude_compiled:
            if regex.search(_scope_text(record, "title_abstract")):
                reasons.append(name)
        for name, regex, scope in require_compiled:
            if not regex.search(_scope_text(record, scope)):
                reasons.append(name)
        types = set((record.get("publication_types") or "").split(";"))
        if EXCLUDED_PUBLICATION_TYPES & types:
            reasons.append("excluded_publication_type")

        for reason in reasons:
            rejected_by_rule.setdefault(reason, set()).add(record["pmid"])
        if reasons:
            excluded.append({**record, "excluded_by": ";".join(reasons)})
        else:
            kept.append(record)
        progress.update()
    progress.close()

    union = set().union(*rejected_by_rule.values()) if rejected_by_rule else set()
    assert {r["pmid"] for r in excluded} == union, \
        "rule application is not the union of the per-rule rejection sets"

    all_rules = ([n for n, _ in EXCLUDE_RULES]
                 + [n for n, _, _ in REQUIRE_RULES]
                 + ["excluded_publication_type"])
    stats = [
        {
            "rule": name,
            "records_rejected": len(rejected_by_rule.get(name, ())),
            "rejected_only_by_this_rule": sum(
                1 for r in excluded if r["excluded_by"] == name),
        }
        for name in all_rules
    ]
    return kept, excluded, stats


def write_csv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_reading_list(path: str, records: list[dict]) -> None:
    ordered = sorted(records, key=lambda r: (r.get("year", ""), r.get("title", "")))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Candidate records for manual screening\n")
        fh.write(f"Generated {datetime.date.today().isoformat()} -- "
                 f"{len(ordered)} records\n")
        fh.write("=" * 79 + "\n\n")
        for index, record in enumerate(ordered, 1):
            fh.write(f"{index:4d}. {record.get('title', '')}\n")
            fh.write(f"      {record.get('journal', '')} ({record.get('year', '')})"
                     f"  PMID {record.get('pmid', '')}")
            if record.get("doi"):
                fh.write(f"  doi:{record['doi']}")
            fh.write("\n\n")


def write_summary(path: str, strand_hits: dict, n_union: int, n_kept: int,
                  n_excluded: int, stats: list[dict]) -> None:
    lines = []
    add = lines.append
    add("Reproducible literature search -- summary")
    add("=" * 79)
    add(f"Run date            : {datetime.date.today().isoformat()}")
    add("Database            : PubMed (NCBI E-utilities)")
    add("Date range in query : 2016-01-01 to 2026-12-31, English")
    add("")
    add("Stage 1 -- search strands")
    add("-" * 79)
    for name, hits in strand_hits.items():
        add(f"  {name:38s} {hits:>7d} records")
    add(f"  {'UNION (deduplicated by PubMed)':38s} {n_union:>7d} records")
    add("")
    add("Stage 3 -- automated screening")
    add("-" * 79)
    add(f"  {'Records entering screening':44s} {n_union:>7d}")
    for stat in stats:
        add(f"    rejected by {stat['rule']:30s} {stat['records_rejected']:>7d}"
            f"   (sole reason: {stat['rejected_only_by_this_rule']})")
    add(f"  {'Records removed (union of all rules)':44s} {n_excluded:>7d}")
    add(f"  {'CANDIDATES FOR MANUAL SCREENING':44s} {n_kept:>7d}")
    add("")
    add("Reference run by the authors")
    add("-" * 79)
    ref = REFERENCE_RUN
    add(f"  Date                              : {ref['date']}")
    add(f"  Union records                     : {ref['union_records']}")
    add(f"  Candidates after screening        : {ref['candidates_after_screening']}")
    add(f"  Studies selected after full text  : {ref['studies_selected_manually']}")
    add("")
    add("PubMed grows continuously, so a later run is expected to return more")
    add("records. The queries and the rules are fixed, so a difference in counts")
    add("reflects database growth rather than a change of method.")
    add("")
    add("Stage 1 -- full query strings")
    add("=" * 79)
    for name, query in STRANDS.items():
        add("")
        add(f"{name}:")
        add(query)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# =============================================================================
# Offline tests of the screening logic (no network required)
# =============================================================================

def selftest() -> None:
    records = [
        {"pmid": "1", "title": "Deep learning of TILs in colorectal cancer",
         "abstract": "We analysed whole-slide images from 500 patients.",
         "publication_types": "Journal Article"},
        {"pmid": "2", "title": "Radiomics of rectal cancer",
         "abstract": "Magnetic resonance imaging features of the immune stroma.",
         "publication_types": "Journal Article"},
        {"pmid": "3", "title": "Gene signature in colon cancer",
         "abstract": "Bulk RNA sequencing of 300 immune-infiltrated tumours.",
         "publication_types": "Journal Article"},
        {"pmid": "4", "title": "A review of immune markers",
         "abstract": "Stained images of colorectal tissue are discussed.",
         "publication_types": "Journal Article;Review"},
        {"pmid": "5", "title": "CD3 and CD8 immunohistochemistry in the QUASAR trial",
         "abstract": ("Pathologists, assisted by artificial intelligence, calculated "
                      "CD3 and CD8 cell densities in colorectal tumour sections."),
         "publication_types": "Journal Article;Randomized Controlled Trial"},
        {"pmid": "6", "title": "Knockdown of a gene in colorectal cells",
         "abstract": ("Western blot and qPCR after siRNA transfection; immune "
                      "staining of sections was also scored."),
         "publication_types": "Journal Article"},
    ]
    kept, excluded, stats = apply_rules(records)
    kept_ids = {r["pmid"] for r in kept}
    reasons = {r["pmid"]: r["excluded_by"] for r in excluded}

    assert kept_ids == {"1", "5"}, f"unexpected kept set: {sorted(kept_ids)}"
    assert "other_imaging_modality" in reasons["2"]
    assert "tissue_observation_required" in reasons["3"]
    assert "excluded_publication_type" in reasons["4"]
    assert "wet_lab_assay" in reasons["6"]
    # a record may fail several rules; every one of them is recorded
    assert ";" in reasons["4"], "multiple rejection reasons must all be listed"
    # record 5 is the QUASAR-style study: no imaging word, kept via
    # "immunohistochemistry" and "sections"
    assert "5" not in reasons

    union = set()
    for record in excluded:
        union.add(record["pmid"])
    assert union == set(reasons), "excluded set is not the union of rejections"

    print("selftest passed")
    print(f"  kept     : {sorted(kept_ids)}")
    for pmid, why in sorted(reasons.items()):
        print(f"  excluded : {pmid}  <- {why}")


# =============================================================================

def main() -> None:
    global _log_path, _quiet

    parser = argparse.ArgumentParser(
        description="Reproduce the candidate corpus for the CRC-TME review.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=None, help="output directory")
    parser.add_argument("--no-cache", action="store_true",
                        help="re-download even if a harvest exists")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress bars")
    parser.add_argument("--selftest", action="store_true",
                        help="run offline tests of the screening logic and exit")
    args = parser.parse_args()

    _quiet = args.quiet
    if args.selftest:
        selftest()
        return

    if not os.environ.get("NCBI_EMAIL"):
        sys.exit('Set NCBI_EMAIL before running, as required by the NCBI usage '
                 'policy:\n    export NCBI_EMAIL="your.address@institution.org"')

    outdir = args.outdir or f"search_{datetime.date.today().isoformat()}"
    os.makedirs(outdir, exist_ok=True)
    _log_path = os.path.join(outdir, "run_log.txt")
    log(f"Output directory: {os.path.abspath(outdir)}")
    if not _HAVE_TQDM and not _quiet:
        log("(install 'tqdm' for nicer progress bars; using the built-in one)")

    log("Stage 1: running the search strands")
    strand_hits = {}
    for name, query in STRANDS.items():
        hits = count_hits(query)
        strand_hits[name] = hits
        log(f"  {name}: {hits} records")

    union_path = os.path.join(outdir, "union_records.csv")
    fields = ["pmid", "year", "journal", "title", "abstract",
              "publication_types", "doi", "record_type"]
    if os.path.exists(union_path) and not args.no_cache:
        records = list(csv.DictReader(open(union_path, encoding="utf-8")))
        log(f"Stage 2: reusing the existing harvest ({len(records)} records); "
            f"pass --no-cache to download again")
    else:
        log("Stage 2: harvesting the union")
        records = harvest_union()
        write_csv(union_path, records, fields)
        log(f"  wrote {union_path} ({len(records)} records)")

    log("Stage 3: applying the screening rules")
    kept, excluded, stats = apply_rules(records)
    for stat in stats:
        log(f"  {stat['rule']}: rejected {stat['records_rejected']} "
            f"(sole reason for {stat['rejected_only_by_this_rule']})")

    write_csv(os.path.join(outdir, "candidate_records.csv"), kept, fields)
    write_csv(os.path.join(outdir, "excluded_records.csv"), excluded,
              fields + ["excluded_by"])
    write_reading_list(os.path.join(outdir, "candidate_records.txt"), kept)
    write_summary(os.path.join(outdir, "search_summary.txt"),
                  strand_hits, len(records), len(kept), len(excluded), stats)

    log("")
    log(f"Union records                 : {len(records)}")
    log(f"Removed by screening rules    : {len(excluded)}")
    log(f"CANDIDATES for manual review  : {len(kept)}")
    ref = REFERENCE_RUN
    if ref["union_records"]:
        log(f"(authors' reference run: {ref['union_records']} -> "
            f"{ref['candidates_after_screening']} on {ref['date']})")
    else:
        log(f"(authors' reference run: {ref['candidates_after_screening']} "
            f"candidates on {ref['date']})")
    log("")
    log(f"Main output: {os.path.join(outdir, 'candidate_records.csv')}")
    log(f"Summary    : {os.path.join(outdir, 'search_summary.txt')}")


if __name__ == "__main__":
    main()
