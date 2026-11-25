"""
Central path configuration for the H9 diploid ASE pipeline.

Every script imports its locations from here, so nothing is tied to one
person's machine. Two environment variables control it:

  H9_ASE_DATA_ROOT   Working directory holding the large inputs and
                     intermediate outputs (genomes/, ase/, ase_v2/,
                     annotations/, results/, vep_redo/, variant_catalog/,
                     clean_analysis/). Default: ~/reganalysis
  SAMTOOLS           samtools executable. Default: the one on PATH.

Example:
  export H9_ASE_DATA_ROOT=/data/h9_ase
  python scripts/paper_published/01_build_master_variants.py
"""

import os
import shutil
import sys

# Repository root (the directory containing scripts/, data/, results/, figures/)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Working data directory (large files that are not committed to the repo)
DATA_ROOT = os.path.abspath(os.path.expanduser(
    os.environ.get('H9_ASE_DATA_ROOT') or '~/reganalysis'))

# External tools
SAMTOOLS = os.environ.get('SAMTOOLS') or shutil.which('samtools') or 'samtools'

if not os.path.isdir(DATA_ROOT):
    print(f"WARNING: data directory not found: {DATA_ROOT}\n"
          f"         Set H9_ASE_DATA_ROOT to the folder holding the pipeline "
          f"inputs (see data/EXTERNAL_DATA.md).", file=sys.stderr)


def data_path(*parts):
    """Path inside the working data directory."""
    return os.path.join(DATA_ROOT, *parts)


def repo_path(*parts):
    """Path inside this repository."""
    return os.path.join(REPO_ROOT, *parts)
