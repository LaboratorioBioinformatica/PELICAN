#!/usr/bin/env python3
"""
pelican-setup: one-time preparation and installation check for PELICAN.

- checks the external tools PELICAN calls
- installs a GeneMark license key (--gm-key)
- tests GeneMarkS on a small reference genome
- downloads/prepares the COG + PHROG BLAST database (same files PELICAN uses at run time)
- downloads the DNABERT-6 model to the HuggingFace cache
"""
import argparse
import os
import shutil
import sys
import tempfile

from pelican.genemarks.podman_genemarks import run_podman as gm

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
DATABASES_DIR = os.path.join(PACKAGE_DIR, "databases")

REQUIRED_TOOLS = {
    "prokka": "conda install -c conda-forge -c bioconda prokka",
    "aragorn": "conda install -c conda-forge -c bioconda aragorn",
    "psiblast": "conda install -c conda-forge -c bioconda blast",
    "makeblastdb": "conda install -c conda-forge -c bioconda blast",
    "blastp": "conda install -c conda-forge -c bioconda blast",
    "blastn": "conda install -c conda-forge -c bioconda blast",
    "wget": "conda install -c conda-forge wget",
    "perl": "conda install -c conda-forge perl",
    "phanotate.py": "pip install phanotate",
}

results = []


def report(name, ok, detail=""):
    results.append((name, ok, detail))
    status = "OK  " if ok else "FAIL"
    print(f"[{status}] {name}" + (f" -- {detail}" if detail else ""))


def check_tools():
    for tool, fix in REQUIRED_TOOLS.items():
        path = shutil.which(tool)
        report(f"tool: {tool}", path is not None, path or f"not found. Install with: {fix}")


def install_gm_key(key_path):
    """Copy a GeneMark key to ~/.gm_key (the standard GeneMark location)."""
    if not os.path.isfile(key_path):
        report("GeneMark key install", False, f"file not found: {key_path}")
        return
    target = os.path.expanduser("~/.gm_key")
    if os.path.exists(target):
        shutil.copy(target, target + ".bak")
    with open(target, "wb") as fh:
        fh.write(gm._read_key(key_path))
    report("GeneMark key install", True, f"installed to {target}")


def check_genemark():
    backend = gm.prepare_genemark()
    if backend is None:
        report("GeneMarkS backend", False, "no usable backend (native, podman or docker)")
        return
    report("GeneMarkS backend", True, backend)

    fasta = os.path.join(gm.SELFTEST_DIR, "NC_010326.1.fasta")
    expected = os.path.join(gm.SELFTEST_DIR, "NC_010326.1.fasta.lst.expected")
    with tempfile.TemporaryDirectory(prefix="pelican_selftest_") as tmp:
        gm.run_podman(fasta, tmp)
        produced = os.path.join(tmp, "NC_010326.1.fasta.lst")
        if gm.license_expired(os.path.join(tmp, gm.LOG_NAME)):
            report("GeneMarkS test run", False, "license key expired (see instructions above)")
            return
        if not os.path.exists(produced):
            log = os.path.join(tmp, gm.LOG_NAME)
            tail = open(log, errors="ignore").read()[-800:] if os.path.exists(log) else ""
            report("GeneMarkS test run", False, f"no prediction produced. Log tail:\n{tail}")
            return

        def genes(path):
            with open(path) as fh:
                return [l for l in fh if not l.startswith("Date:")]
        same = genes(produced) == genes(expected)
        report("GeneMarkS test run", same,
               "predictions identical to reference" if same else "predictions differ from reference output")


def prepare_databases():
    if not os.access(DATABASES_DIR, os.W_OK):
        report("database folder writable", False,
               f"{DATABASES_DIR} is not writable. Install PELICAN in a user-owned environment (e.g. conda).")
        return
    # Same helpers and paths used by PELICAN at run time (pelican/core.py)
    from pelican import core
    phrogs_db = os.path.join(DATABASES_DIR, "Phrogs_genes.faa")
    cogs_db = os.path.join(DATABASES_DIR, "COGorg24.faa")
    combined_db = os.path.join(DATABASES_DIR, "combined_db.faa")
    try:
        core.download_COG(cogs_db)
        core.ensure_database_uncompressed(cogs_db)
        core.ensure_database_uncompressed(phrogs_db)
        if not os.path.exists(combined_db):
            core.concat_fasta_files([phrogs_db, cogs_db], combined_db)
        if not os.path.exists(f"{combined_db}.pin"):
            core.make_blast_db(combined_db)
    except Exception as e:
        report("BLAST database", False, str(e))
        return
    ok = os.path.exists(f"{combined_db}.pin")
    report("BLAST database", ok, combined_db if ok else "makeblastdb did not create the database")


def prepare_dnabert():
    from pelican import core
    try:
        core.load_dnabert6()
        report("DNABERT-6 model", True, f"cached ({core.MODEL_ID})")
    except Exception as e:
        report("DNABERT-6 model", False, f"{e}. Check internet access or set HF_TOKEN / HF_HOME.")


def main():
    parser = argparse.ArgumentParser(description="Prepare and check a PELICAN installation")
    parser.add_argument("--gm-key", help="GeneMark license key file (gm_key_64.gz or plain) to install as ~/.gm_key")
    parser.add_argument("--check-only", action="store_true", help="Only check tools and GeneMarkS; do not download anything")
    parser.add_argument("--skip-db", action="store_true", help="Do not prepare the BLAST database")
    parser.add_argument("--skip-dnabert", action="store_true", help="Do not download the DNABERT-6 model")
    args = parser.parse_args()

    print("=" * 60)
    print("PELICAN setup")
    print("=" * 60)
    if args.gm_key:
        install_gm_key(args.gm_key)
    check_tools()
    check_genemark()
    if not args.check_only:
        if not args.skip_db:
            prepare_databases()
        if not args.skip_dnabert:
            prepare_dnabert()

    failed = [r for r in results if not r[1]]
    print("=" * 60)
    if failed:
        print(f"{len(failed)} problem(s) found:")
        for name, _, detail in failed:
            print(f"  - {name}: {detail.splitlines()[0] if detail else ''}")
        sys.exit(1)
    print("Everything is ready. Try: pelican -i genome.fasta --output_path results")


if __name__ == "__main__":
    main()
