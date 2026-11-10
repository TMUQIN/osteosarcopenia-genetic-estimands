#!/usr/bin/env python3
"""Repository hygiene smoke test only; does not validate scientific estimates."""
from pathlib import Path
import re
ROOT = Path(__file__).resolve().parents[1]
required = [
    "README.md", "CITATION.cff", "CODE_AVAILABILITY.md",
    "DATA_AVAILABILITY.md", "REPRODUCIBILITY_LIMITATIONS.md",
    "code/verification/model_FROZEN.lav",
    "code/verification/model_DERIVED.lav",
    "code/verification/verify_v7_independent.R",
    "data/derived_public_candidate/FIGURE2_DATA.tsv",
    "data/derived_public_candidate/BONE_OUTCOME_SHIFTS_DATA.tsv",
    "data/derived_public_candidate/FIGURE_S1_DATA.tsv",
]
missing = [x for x in required if not (ROOT / x).is_file()]
assert not missing, f"Missing required repository files: {missing}"
for forbidden in ["verify_s9.py", "TABLE_S9", "whole_body_fat_mass.sumstats", "standing_height"]:
    hits = [p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file() and forbidden.lower() in p.name.lower()]
    assert not hits, f"Forbidden/currently excluded material found: {forbidden}: {hits}"
absolute_path_pattern = re.compile(r"(?:[A-Za-z]:[\\/]|/mnt/data/|/home/[^ ]+)")
scan_ext = {".py", ".R", ".lav", ".md", ".txt", ".cff"}
scanned = 0
for p in ROOT.rglob("*"):
    if p.is_file() and p != Path(__file__) and p.suffix in scan_ext:
        scanned += 1
        text = p.read_text(encoding="utf-8", errors="ignore")
        assert not absolute_path_pattern.search(text), f"Possible local absolute path in {p.relative_to(ROOT)}"

# Validate the release inventory against the current working tree (excluding itself).
import csv, hashlib
inventory_path = ROOT / "PUBLIC_FILE_INVENTORY.tsv"
with inventory_path.open(encoding="utf-8", newline="") as f:
    inventory_rows = list(csv.DictReader(f, delimiter="\t"))
listed = {row["relative_path"]: row for row in inventory_rows}
actual = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob("*") if p.is_file() and p != inventory_path}
assert set(listed) == set(actual), f"Inventory path mismatch: missing={sorted(set(actual)-set(listed))}; absent={sorted(set(listed)-set(actual))}"
for rel, path in actual.items():
    payload = path.read_bytes()
    row = listed[rel]
    assert int(row["bytes"]) == len(payload), f"Inventory size mismatch: {rel}"
    assert row["sha256"] == hashlib.sha256(payload).hexdigest(), f"Inventory SHA-256 mismatch: {rel}"

print("REPOSITORY_HYGIENE_SMOKE_TEST=PASS")
print(f"required_files={len(required)}")
print(f"text_files_scanned={scanned}")
print("NOTE=This test does not run LDSC, GenomicSEM, or verify scientific estimates.")
