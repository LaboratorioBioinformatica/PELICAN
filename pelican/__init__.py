"""
PELICAN - Phage gEnome tooL for Inference of Consensus ANotation

A Machine Learning based tool for automated genome consensus annotation.
"""

__version__ = "0.1.5"
__author__ = "Fernando Rossi"
__email__ = "fenobrerossi@gmail.com"

from .core import run_pelican, main

__all__ = ["run_pelican", "main"]
