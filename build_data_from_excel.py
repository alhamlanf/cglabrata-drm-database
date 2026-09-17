#!/usr/bin/env python3
# =============================================================================
# build_data_from_excel.py
# -----------------------------------------------------------------------------
# Regenerates the four CSVs in data/ from the master Excel workbook.
# Run:  python3 build_data_from_excel.py
#
# Design decisions (locked with the user, see PROJECT_BRIEF.md / README.md):
#   * Mutations are rebuilt keyed by Excel ROW NUMBER:  mutation_id = MUT_R{row}.
#     (Rows may be inserted/reordered between Excel versions; IDs follow rows.)
#   * The per-reference discordance phenotype table was hand-curated and cannot
#     be regenerated from the Excel. So for a mutation whose NAME still matches a
#     mutation in the previous CSV, we CARRY FORWARD its curated per-reference
#     R/I/S + detail/comment. For new or renamed mutations, every cited reference
#     gets the row-level AGGREGATE phenotype (Excel cols "ANF/CAS/MCF <drug>").
#   * Reference metadata for IDs not already in references.csv is pulled from
#     References_list_copy_v2.xlsx; existing reference rows are kept verbatim.
#   * The 4 legacy aliases are re-keyed to the new IDs by the target molecule
#     named in each alias note (the old alias IDs never matched by row).
#
# openpyxl is used (not R/readxl) because it is the only convenient way to read
# the workbook's threaded comments and cached values.
# =============================================================================

import csv
import os
import re
import shutil
import sys
import warnings
from collections import Counter, defaultdict
from datetime import datetime

import openpyxl

warnings.filterwarnings("ignore")  # silence "Data Validation extension" notice

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
EXCEL = os.path.join(HERE, "Database for C.glabrata DRM_clean5.xlsx")
REF_META_XLSX = os.path.join(os.path.dirname(HERE), "References_list_copy_v2.xlsx")
MAIN_SHEET = "DRMs_Echinocandins"

# Column positions in the main sheet (1-indexed).
C_GENE, C_HOTSPOT, C_COMB, C_NAME, C_TYPE, C_VARPHENO = 1, 2, 3, 4, 5, 6
C_ANF_AGG, C_CAS_AGG, C_MCF_AGG = 10, 11, 12   # row-level aggregate phenotype
C_REF_START = 13                                # reference-id list begins here

# Aggregate-phenotype word -> R/I/S code.
AGG_MAP = {
    "resistant": "R",
    "intermediate": "I",
    "sensitive": "S",
    "no phenotypic test": "",
    "no data": "",
    "": "",
}

# Legacy aliases re-keyed to the new row-based IDs by target molecule.
# (old alias id -> new mutation_id for the molecule the alias note refers to)
ALIAS_RELINK = {
    "S659Y":               "MUT_R050",  # F659Y (pos 659 is Phe, not Ser)
    "I660_L661insF":       "MUT_R241",  # clean5 L660_I661insF (letters inverted)
    ("S645P (C. albicans numbering)", "FKS1"): "MUT_R008",  # S629P homolog
    ("S645P (C. albicans numbering)", "FKS2"): "MUT_R003",  # S663P homolog
}


# ---- helpers ----------------------------------------------------------------
def norm(v):
    """Trim, drop non-breaking spaces."""
    return "" if v is None else str(v).replace("\xa0", " ").strip()


def name_key(name):
    """Canonical form for matching a mutation across Excel versions."""
    return norm(name).upper().replace("STOP", "*")


def parse_residue(name):
    """Return (wt, position, mutant) parsed from a single mutation name.
    Blank triple for anything that doesn't look like a single substitution."""
    m = re.match(r"^([A-Za-z])(\d+)(.*)$", norm(name))
    if not m:
        return "", "", ""
    wt, pos, rest = m.group(1), m.group(2), m.group(3).strip()
    return wt, pos, rest


def read_existing_csv(fname):
    path = os.path.join(DATA, fname)
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# =============================================================================
# 1. Load the previous CSVs (for carry-forward + reference reuse)
# =============================================================================
old_mutations = read_existing_csv("mutations.csv")
old_pheno_rows = read_existing_csv("mutation_reference_phenotype.csv")
old_refs = read_existing_csv("references.csv")

# name_key -> old mutation_id, preferring the id that has the most phenotype rows.
pheno_count = Counter(r["mutation_id"] for r in old_pheno_rows)
old_id_by_name = {}
best_count = {}
for r in old_mutations:
    k = name_key(r["mutation_name"])
    c = pheno_count.get(r["mutation_id"], 0)
    if k not in old_id_by_name or c > best_count[k]:
        old_id_by_name[k] = r["mutation_id"]
        best_count[k] = c

# old_id -> {reference_id -> full phenotype row}
old_pheno = defaultdict(dict)
for r in old_pheno_rows:
    old_pheno[r["mutation_id"]][r["reference_id"]] = r

# =============================================================================
# 2. Read the master Excel main sheet
# =============================================================================
if not os.path.exists(EXCEL):
    sys.exit("Cannot find workbook: " + EXCEL)
wb = openpyxl.load_workbook(EXCEL, data_only=True)
ws = wb[MAIN_SHEET]

last_row = 1
for r in range(2, ws.max_row + 1):
    if norm(ws.cell(r, C_NAME).value):
        last_row = r

# First pass: collect single-mutation names -> new id (for combo components).
new_id_by_name = {}
for r in range(2, last_row + 1):
    nm = norm(ws.cell(r, C_NAME).value)
    if nm and norm(ws.cell(r, C_COMB).value) != "Y" and "/" not in nm:
        new_id_by_name.setdefault(name_key(nm), "MUT_R%03d" % r)

# =============================================================================
# 3. Build mutations.csv  and  mutation_reference_phenotype.csv
# =============================================================================
mutations_out = []
pheno_out = []
cited_refs = set()
unknown_agg = Counter()
carried, aggregated = 0, 0

for r in range(2, last_row + 1):
    name = norm(ws.cell(r, C_NAME).value)
    if not name:
        continue
    mid = "MUT_R%03d" % r
    gene = norm(ws.cell(r, C_GENE).value)
    hotspot = norm(ws.cell(r, C_HOTSPOT).value)
    comb = norm(ws.cell(r, C_COMB).value)
    mtype = norm(ws.cell(r, C_TYPE).value)
    varpheno = norm(ws.cell(r, C_VARPHENO).value)

    is_comb = "Y" if (comb == "Y" or (comb == "" and "/" in name)) else "N"

    # Combination mutations span multiple residues (and often both genes), so a
    # single gene/hot-spot label doesn't apply — they are labelled "multiple".
    # Enforced here so the rule holds even if an Excel row is missed.
    if is_comb == "Y":
        gene = "multiple"
        hotspot = "multiple"

    # Residue fields only for single substitutions.
    if is_comb == "N":
        wt, pos, mut = parse_residue(name)
    else:
        wt, pos, mut = "", "", ""

    # Combination component ids (best effort; not used by the app).
    comp = ""
    if is_comb == "Y":
        parts = [p.strip() for p in name.split("/") if p.strip()]
        comp = ";".join(new_id_by_name.get(name_key(p), "?" + p) for p in parts)

    mutations_out.append({
        "mutation_id": mid, "excel_row": r, "gene": gene, "hotspot": hotspot,
        "mutation_name": name, "is_combination": is_comb, "wt_residue": wt,
        "position": pos, "mutant_residue": mut, "mutation_type": mtype,
        "variable_phenotype": varpheno, "component_mutation_ids": comp,
    })

    # ---- reference list for this mutation ----
    reflist = []
    for c in range(C_REF_START, ws.max_column + 1):
        v = ws.cell(r, c).value
        if isinstance(v, int):
            reflist.append(str(v))
        elif norm(v).isdigit():
            reflist.append(norm(v))
    cited_refs.update(reflist)

    # ---- aggregate phenotype for this row ----
    def agg(col):
        raw = norm(ws.cell(r, col).value).lower()
        if raw not in AGG_MAP:
            unknown_agg[raw] += 1
        return AGG_MAP.get(raw, "")
    a_anf, a_cas, a_mcf = agg(C_ANF_AGG), agg(C_CAS_AGG), agg(C_MCF_AGG)

    oldid = old_id_by_name.get(name_key(name))
    for ref in reflist:
        if oldid and ref in old_pheno.get(oldid, {}):
            # carry forward the hand-curated per-reference record
            src = old_pheno[oldid][ref]
            pheno_out.append({
                "mutation_id": mid, "reference_id": ref,
                "anf_class": src["anf_class"], "cas_class": src["cas_class"],
                "mcf_class": src["mcf_class"], "anf_detail": src["anf_detail"],
                "cas_detail": src["cas_detail"], "mcf_detail": src["mcf_detail"],
                "anf_comment": src["anf_comment"], "cas_comment": src["cas_comment"],
                "mcf_comment": src["mcf_comment"],
            })
            carried += 1
        else:
            # fall back to the row-level aggregate
            pheno_out.append({
                "mutation_id": mid, "reference_id": ref,
                "anf_class": a_anf, "cas_class": a_cas, "mcf_class": a_mcf,
                "anf_detail": "", "cas_detail": "", "mcf_detail": "",
                "anf_comment": "", "cas_comment": "", "mcf_comment": "",
            })
            aggregated += 1

# =============================================================================
# 4. Build references.csv  (keep existing rows, add missing cited ids)
# =============================================================================
ref_by_id = {r["reference_id"]: r for r in old_refs}
missing = sorted((int(x) for x in cited_refs if x not in ref_by_id))
ref_cols = ["reference_id", "pmid", "doi", "pmc", "first_author", "year",
            "journal", "title", "study_type", "full_citation", "authors"]

if missing:
    wbref = openpyxl.load_workbook(REF_META_XLSX, data_only=True)
    wsref = wbref[wbref.sheetnames[0]]
    hdr = {norm(c.value): i + 1 for i, c in enumerate(wsref[1])}
    meta = {}
    for rr in range(2, wsref.max_row + 1):
        idv = wsref.cell(rr, hdr["#"]).value
        if isinstance(idv, int):
            g = lambda key: norm(wsref.cell(rr, hdr[key]).value) if key in hdr else ""
            meta[idv] = {
                "reference_id": str(idv), "pmid": g("PMID"), "doi": g("DOI"),
                "pmc": "", "first_author": g("First Author"),
                "year": g("Publication Year"), "journal": g("Journal/Book"),
                "title": g("Title"), "study_type": "",
                "full_citation": g("Citation"), "authors": g("Authors"),
            }
    for mid_ in missing:
        if mid_ in meta:
            ref_by_id[str(mid_)] = meta[mid_]
        else:
            ref_by_id[str(mid_)] = {c: "" for c in ref_cols}
            ref_by_id[str(mid_)]["reference_id"] = str(mid_)

refs_out = [ref_by_id[k] for k in sorted(ref_by_id, key=lambda x: int(x))
            if k in {str(int(x)) for x in cited_refs}]

# =============================================================================
# 5. Build mutation_aliases.csv  (re-keyed to new IDs)
# =============================================================================
old_aliases = read_existing_csv("mutation_aliases.csv")
aliases_out = []
for a in old_aliases:
    label = a["alias_label"]
    note = a["notes"]
    if label in ALIAS_RELINK:
        newid = ALIAS_RELINK[label]
    else:
        # the two "S645P (C. albicans numbering)" aliases differ by target gene,
        # identified from the homolog molecule named in the note
        # (FKS1 homolog -> S629P, FKS2 homolog -> S663P). Note text mentions
        # "C. albicans FKS1 numbering" for both, so key on the homolog name.
        gene = "FKS2" if "S663P" in note else "FKS1"
        newid = ALIAS_RELINK.get((label, gene), a["mutation_id"])
    aliases_out.append({"mutation_id": newid, "alias_label": label, "notes": note})

# =============================================================================
# 6. Back up old CSVs and write new ones
# =============================================================================
backup = os.path.join(HERE, "data_backups",
                      "backup_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
os.makedirs(backup, exist_ok=True)
for fn in ["mutations.csv", "mutation_reference_phenotype.csv",
           "references.csv", "mutation_aliases.csv"]:
    shutil.copy2(os.path.join(DATA, fn), os.path.join(backup, fn))


def write_csv(fname, rows, cols):
    with open(os.path.join(DATA, fname), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


write_csv("mutations.csv", mutations_out,
          ["mutation_id", "excel_row", "gene", "hotspot", "mutation_name",
           "is_combination", "wt_residue", "position", "mutant_residue",
           "mutation_type", "variable_phenotype", "component_mutation_ids"])
write_csv("mutation_reference_phenotype.csv", pheno_out,
          ["mutation_id", "reference_id", "anf_class", "cas_class", "mcf_class",
           "anf_detail", "cas_detail", "mcf_detail",
           "anf_comment", "cas_comment", "mcf_comment"])
write_csv("references.csv", refs_out, ref_cols)
write_csv("mutation_aliases.csv", aliases_out,
          ["mutation_id", "alias_label", "notes"])

# =============================================================================
# 7. Report
# =============================================================================
print("Backup of previous CSVs -> %s" % backup)
print("mutations.csv                    : %d rows" % len(mutations_out))
print("mutation_reference_phenotype.csv : %d rows (%d carried-forward, %d aggregate)"
      % (len(pheno_out), carried, aggregated))
print("references.csv                   : %d rows (added %s)"
      % (len(refs_out), missing or "none"))
print("mutation_aliases.csv             : %d rows -> %s"
      % (len(aliases_out), [a["mutation_id"] for a in aliases_out]))
if unknown_agg:
    print("NOTE unrecognised aggregate phenotype words:", dict(unknown_agg))
