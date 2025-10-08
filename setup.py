from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements from the dependencies
requirements = [
    "pandas>=1.3.0",
    "matplotlib>=3.3.0",
    "pyrodigal>=2.0.0",
    "biopython>=1.79",
    "tqdm>=4.62.0",
    "scikit-learn>=1.0.0",
    "mycolorpy>=1.5.0",
    "numpy>=1.21.0",
    "scipy>=1.7.0",
    "pygenomeviz>=0.4.0",
    "upsetplot>=0.6.0",
    "bcbio-gff>=0.6.9",
    "protpy>=1.0.0",
    "joblib>=1.1.0",
    "torch>=1.9.0",
    "h5py>=3.3.0",
    "transformers>=4.12.0",
]

setup(
    name="pelican-annotation",
    version="0.1.5",
    author="Fernando Rossi",
    author_email="fenobrerossi@gmail.com",
    description="Phage gEnome tooL for Inference of Consensus ANotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FefoRossi/PELICAN",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "pelican=pelican.cli:cli_main",
        ],
    },
    include_package_data=True,
    package_data={
        "pelican": [
            "databases/*",
            "models/*",
            "genemarks/*",
        ],
    },
    zip_safe=False,
)
