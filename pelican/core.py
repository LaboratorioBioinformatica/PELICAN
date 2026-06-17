def concat_fasta_files(input_files, output_file):
# Concatenate multiple FASTA files into one output file.
    """Concatena múltiplos arquivos FASTA em um único arquivo de saída."""
    with open(output_file, 'w') as outfile:
        for fname in input_files:
            with open(fname, 'r') as infile:
                shutil.copyfileobj(infile, outfile)

def make_blast_db(fasta_file, db_type='prot'):
# Create a BLAST database from a FASTA file.
    print(color_text(f"Creating BLAST database", "yellow"))
    cmd = f"makeblastdb -in {fasta_file} -dbtype {db_type} -logfile /dev/null"
    subprocess.run(cmd, shell=True, check=True)

def codon_usage(seq):
# Calculate codon usage frequency for a given sequence.
    seq = seq.upper().replace('T', 'U')
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    total = len(codons)
    counts = {codon: 0 for codon in valid_codons}
    for codon in codons:
        if codon in counts:
            counts[codon] += 1
    # Frequência relativa
    freq_vec = [counts[codon]/total if total > 0 else 0 for codon in valid_codons]
    # DataFrame: códons como colunas, valores como frequência relativa
    df_cod = pd.DataFrame([freq_vec], columns=valid_codons)
    return df_cod

def parse_predictions(output_file: str) -> list:
# Parse gene prediction results from GeneMarkS-2 output.
        """
        Faz parsing dos resultados de predição no formato LST do GeneMarkS-2
        Args:
            output_file: Arquivo de saída do GeneMarkS-2
        Returns:
            Lista de dicionários com informações dos genes preditos
        """
        import re
        genes = []
        current_seqid = None
        fallback_genes = []
        splitter = re.compile(r"\s+")
        num_cleaner = r"[^0-9]"
        with open(output_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('SequenceID:'):
                    current_seqid = line.split(':', 1)[1].strip()
                    continue
                parts = splitter.split(line)
                if not parts or not parts[0].isdigit():
                    continue
                if len(parts) >= 9:
                    try:
                        gene_id = int(parts[0])
                        if parts[1] not in ['+', '-']:
                            continue
                        # Interpretar corretamente >N e <N
                        def parse_pos(val):
                            if val.startswith('>') and val[1:].isdigit():
                                return int(val[1:])
                            elif val.startswith('<') and val[1:].isdigit():
                                return int(val[1:])
                            else:
                                digits = re.sub(r'[^0-9]', '', val)
                                return int(digits) if digits else 0
                        start = parse_pos(parts[2])
                        end = parse_pos(parts[3])
                        length = int(parts[4])
                        gene = {
                            'gene_id': gene_id,
                            'strand': parts[1],
                            'start': start,
                            'end': end,
                            'length': length,
                            'type': parts[5],
                            'motif': parts[6],
                            'motif_score': int(parts[7]) if parts[7].isdigit() else 0,
                            'frame': int(parts[8]) if parts[8].isdigit() else 1,
                            'seqid': current_seqid
                        }
                        genes.append(gene)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Skipping line {line_num} due to parsing error: {line}")
                        print(f"Error details: {e}")
                        continue
                elif len(parts) == 6:
                    try:
                        gene_id = int(parts[0])
                        strand = parts[1]
                        # Interpretar corretamente >N e <N
                        def parse_pos(val):
                            if val.startswith('>') and val[1:].isdigit():
                                return int(val[1:])
                            elif val.startswith('<') and val[1:].isdigit():
                                return int(val[1:])
                            else:
                                digits = re.sub(r'[^0-9]', '', val)
                                return int(digits) if digits else 0
                        start = parse_pos(parts[2])
                        end = parse_pos(parts[3])
                        length = int(parts[4])
                        gene_class = parts[5]
                        fallback_genes.append({
                            'gene_id': gene_id,
                            'strand': strand,
                            'start': start,
                            'end': end,
                            'length': length,
                            'class': gene_class,
                            'seqid': current_seqid
                        })
                    except Exception as e:
                        print(f"Fallback parse error at line {line_num}: {line}")
                        print(f"Error details: {e}")
                        continue
        if genes:
            return genes
        elif fallback_genes:
            return fallback_genes
        else:
            return []
#!/usr/bin/env python3
"""
PELICAN - Phage Gene Prediction and Annotation Pipeline

COORDINATE SYSTEM NOTES:
- All biological coordinates (from GFF, Prodigal, Phanotate, etc.) are 1-based
- Python string indexing is 0-based
- When extracting sequences from genome string, always convert 1-based to 0-based
- Use helper functions coordinate_to_python_index() and extract_sequence_from_genome()
  to ensure consistent coordinate handling
"""

import pandas as pd
import os
import json
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sys
# sys.stderr = open(os.devnull, 'w')
sys.path.append(os.path.join(os.path.dirname(__file__), '../genemarks'))

from unittest.mock import MagicMock
import types

from pelican.genemarks.python_wrapper import GeneMarkS2Wrapper
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# import seaborn as sns
import pyrodigal_gv
from Bio import SeqIO
import subprocess
from tqdm import tqdm
import io
from Bio.Seq import Seq
import textwrap
import argparse
import random
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.metrics.pairwise import cosine_similarity
from mycolorpy import colorlist as mcp
import numpy as np
from scipy.stats import zscore
import gzip
import shutil
import urllib.request
import tarfile
import tempfile
from scipy.spatial.distance import euclidean

## TESTE COM GENEMARKS ##
from pelican.genemarks.podman_genemarks.run_podman import run_podman

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord
import upsetplot
import re

from BCBio import GFF
from operator import itemgetter
import itertools
# import protpy as protpy

import joblib
from pathlib import Path
# Force PyTorch to use CPU only
# import os
os.environ['TRITON_DISABLE'] = '1'  # Disable Triton JIT compiler
# os.environ['CUDA_VISIBLE_DEVICES'] = ''  # Hide CUDA devices
# os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # Force synchronous CUDA operations
os.environ['CUDA_VISIBLE_DEVICES'] = ''

import torch

 # Use single thread for CPU computations
torch.cuda.is_available = lambda: False
torch.cuda.device_count = lambda: 0
torch.cuda.get_device_capability = lambda *args: (0, 0)
torch.is_cuda = False

property_is_cuda = property(lambda self: False)
torch.Tensor.is_cuda = property_is_cuda

import h5py
from transformers import AutoConfig, T5EncoderModel, T5Tokenizer
import torch.nn as nn

from transformers import BertTokenizer, BertModel, AutoTokenizer, AutoModel, BertConfig

# Disable CUDA after importing torch
# if torch.cuda.is_available():
#     torch.cuda.is_available = lambda: False

import warnings
warnings.filterwarnings('ignore')


######################
## Helper functions ##
######################

def ensure_start_less_than_end(start, end, strand):
# Ensure start coordinate is less than end for consistent handling.
    """
    Ensure that start coordinate is always less than end coordinate,
    regardless of strand direction. This is important for consistent 
    coordinate handling throughout the pipeline.
    
    Args:
        start (int): Start coordinate
        end (int): End coordinate  
        strand (int): Strand direction (1 or -1)
        
    Returns:
        tuple: (corrected_start, corrected_end, strand)
    """
    if start > end:
        # Swap coordinates if start > end
        return end, start, strand
    else:
        return start, end, strand

def coordinate_to_python_index(start_1based, end_1based):
# Convert 1-based coordinates to 0-based Python indices.
    """
    Convert 1-based biological coordinates to 0-based Python indexing.
    
    Args:
        start_1based (int): Start position in 1-based coordinates
        end_1based (int): End position in 1-based coordinates (inclusive)
        
    Returns:
        tuple: (start_0based, end_0based) for Python slicing [start:end]
    """
    return start_1based - 1, end_1based

def extract_sequence_from_genome(genome, start_1based, end_1based, strand=1):
# Extract sequence from genome using coordinates and strand.
    """
    Extract sequence from genome using 1-based coordinates.
    
    Args:
        genome (str): Full genome sequence
        start_1based (int): Start position (1-based, inclusive)
        end_1based (int): End position (1-based, inclusive)
        strand (int): 1 for forward, -1 for reverse
        
    Returns:
        str: Extracted sequence (reverse complemented if strand=-1)
    """
    start_0based, end_0based = coordinate_to_python_index(start_1based, end_1based)
    sequence = genome[start_0based:end_0based]
    
    if strand == -1:
        sequence = str(Seq(sequence).reverse_complement())
    
    return sequence

def ensure_database_uncompressed(db_path):
# Ensure database file is uncompressed, decompress if needed.
    """Ensure the database file is uncompressed, decompress if needed."""
    compressed_path = f"{db_path}.gz"
    
    # If the .gz file exists but the uncompressed doesn't, decompress
    if os.path.exists(compressed_path) and not os.path.exists(db_path):
        print(f"Decompressing {compressed_path}...")
        with gzip.open(compressed_path, 'rb') as f_in:
            with open(db_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Database decompressed to {db_path}")
    
    return db_path

def download_COG(db_path,ftp_link="https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/COGorg24.faa.gz"):
# Download COG database if not present.
    cmd = f"wget -P {os.path.dirname(db_path)} {ftp_link}"
    if not os.path.exists(f"{db_path}.gz"):
        print(color_text(f"Downloading COG database...", "yellow"))
        subprocess.run(cmd, shell=True)


def predict_trna_aragorn(input_fasta, output_file):
# Run Aragorn to predict tRNAs in a FASTA file.
    """
    Run Aragorn to predict tRNAs in the input fasta file.
    
    Args:
        input_fasta (str): Path to input fasta file
        output_file (str): Path to output file for Aragorn results
    """
    command = f"aragorn -t -fon -i {input_fasta} -o {output_file}"
    subprocess.run(command, shell=True)

tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}

def color_text(text, color='green'):
# Colorize text output for terminal display.

    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text

def wprint(text):
# Print wrapped text to terminal.

    print(textwrap.fill(text, width=80, initial_indent="\n  ", 
          subsequent_indent="    ", break_on_hyphens=False))

def check_file_exists(file):
# Check if a file exists, exit if not found.

    if not os.path.exists(file):
        wprint(color_text("The specified input file '" + str(file) + "' does not seem to exist :(", "red"))
        print("\n  Exiting for now.\n")
        exit(1)

def min_max(column):
# Normalize a pandas column by its maximum value.
    max_val = column.max()
    return column / max_val

def convert_gff_to_gbk(input_gff, fasta_file, consensus_df):
# Convert GFF annotation to GenBank format using consensus data.

    gbk_file = "%s.gbk" % os.path.splitext(input_gff)[0]

    with open(gbk_file, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        for record in GFF.parse(input_gff, fasta_handler):
            # sequence in each contig (record)
            record.id = str(record.id)
            subset_seqs_df = consensus_df.copy()
            # get all the seqs in the contigs - and drop the index to reset for 0 indexed loop
            subset_seqs = subset_seqs_df["Sequence"].reset_index(drop=True)
            subset_product = subset_seqs_df["Annotation"].reset_index(drop=True)
            # start the loop
            i = 0

            # instantiate record
            record.annotations["molecule_type"] = "DNA"
            record.annotations["date"] = datetime.today()
            record.annotations["topology"] = "linear"
            record.annotations[
                "data_file_division"
            ] = "USP"  # https://github.com/RyanCook94/inphared/issues/22
            # add features to the record
            for feature in record.features:
                # add translation only if CDS, tRNAs don't need translation
                if feature.type == "CDS":
                    # print(feature)
                    # aa = prot_records[i].seq
                    if feature.location.strand == 1:
                        tmp_seq = Seq(subset_seqs[i]).translate(stop_symbol="")
                        feature.qualifiers.update(
                            {"ID": feature.qualifiers["ID"][0],"protein_id": feature.qualifiers["ID"][0], "locus_tag": feature.qualifiers["ID"][0],"product": subset_product[i],"translation": tmp_seq}  # from the aa seq
                        )
                    else:  # reverse strand -1 needs reverse compliment
                        tmp_seq = str(Seq( subset_seqs[i]).translate(stop_symbol=""))
                        feature.qualifiers.update(
                            {"ID": feature.qualifiers["ID"][0],"protein_id": feature.qualifiers["ID"][0], "locus_tag": feature.qualifiers["ID"][0],"product": subset_product[i],"translation": tmp_seq}  # from the aa seq
                        )
                    i += 1
                elif feature.type == "tRNA":
                    # Para tRNAs, apenas adicionar informações básicas sem tradução
                    feature.qualifiers.update(
                        {"ID": feature.qualifiers["ID"][0], "locus_tag": feature.qualifiers["ID"][0], "product": subset_product[i]}
                    )
                    i += 1
            SeqIO.write(record, gbk_handler, "genbank")

def blast_score(ident, cov, evalue):
# Calculate BLAST score from identity, coverage, and e-value.
    score = ((ident/100) + (cov/100))/2
    return score

def autcnn_score(observed_error):
# Calculate AutoCNN score from observed error.
    score = 1 - observed_error
    return score

# Definir a classe CNNAutoencoder igual ao treinamento
class CNNAutoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim=128):  # ✅ Espaço latente REDUZIDO
        super(CNNAutoencoder, self).__init__()
        self.input_dim = input_dim
        self.encoder_conv = nn.Sequential(
            nn.Conv1d(1, 16, kernel_size=5, stride=2, padding=2),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Dropout(0.05),  # ✅ AUMENTADO de 0.2 para 0.4
            
            nn.Conv1d(16, 32, kernel_size=5, stride=2, padding=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Dropout(0.05),  # ✅ AUMENTADO de 0.1 para 0.3
            
            # ✅ REMOVER terceira camada convolucional para simplificar
        )
        
        with torch.no_grad():
            dummy = torch.zeros(1, 1, input_dim)
            conv_out = self.encoder_conv(dummy)
            self.conv_out_shape = conv_out.shape
            conv_flatten_dim = conv_out.numel() // conv_out.shape[0]
            
        self.encoder_fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(conv_flatten_dim, latent_dim),
            nn.ReLU(),
            nn.Dropout(0.05)  # ✅ Dropout no espaço latente
        )
        
        self.decoder_fc = nn.Sequential(
            nn.Linear(latent_dim, conv_flatten_dim),
            nn.ReLU()
        )
        
        self.decoder_deconv = nn.Sequential(
            nn.Unflatten(1, self.conv_out_shape[1:]),
            nn.ConvTranspose1d(32, 16, kernel_size=5, stride=2, padding=2, output_padding=1),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.ConvTranspose1d(16, 1, kernel_size=5, stride=2, padding=2, output_padding=1),
            nn.Sigmoid()  # ✅ RETORNAR Sigmoid para limitar saída [0,1]
        )

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.encoder_conv(x)
        x = self.encoder_fc(x)
        x = self.decoder_fc(x)
        x = self.decoder_deconv(x)
        x = x.squeeze(1)
        
        if x.shape[1] > self.input_dim:
            x = x[:, :self.input_dim]
        elif x.shape[1] < self.input_dim:
            pad = self.input_dim - x.shape[1]
            x = torch.nn.functional.pad(x, (0, pad))
        return x


##########
## RSCU ##
##########

## dictionary that maps codons to amino acids

# Dicionário que mapeia códons para aminoácidos
codon_aa = {
    "UUU":"Phe", "UUC":"Phe",         
    "UCU":"Ser4", "UCC":"Ser4", "UCA":"Ser4", "UCG":"Ser4",
    "AGU":"Ser2", "AGC":"Ser2",
    "CUU":"Leu4", "CUC":"Leu4", "CUA":"Leu4", "CUG":"Leu4",
    "UUA":"Leu2", "UUG":"Leu2",
    
    "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
    "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
    "CGU":"Arg4", "CGC":"Arg4", "CGA":"Arg4", "CGG":"Arg4",
    "AGA":"Arg2", "AGG":"Arg2",
    "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    
    "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
    "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
   
    "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
    "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}

# Lista de códons válidos (excluindo Met, Trp e stops)
valid_codons = [codon for codon in codon_aa if codon not in ['AUG', 'UGG', 'UAA', 'UAG', 'UGA']]

def codon_features(seq):
# Extract codon features, RSCU, and relative weights from sequence.
    seq = seq.upper().replace('T', 'U')
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    total = len(codons)
    counts = {codon: 0 for codon in valid_codons}
    for codon in codons:
        if codon in counts:
            counts[codon] += 1
    # Frequência relativa
    freq_vec = [counts[codon]/total if total > 0 else 0 for codon in valid_codons]
    # DataFrame para RSCU
    df_cod = pd.DataFrame({'Codon': list(counts.keys()), 'Obs_Freq': list(counts.values())})
    df_cod['Amino_Acid'] = [codon_aa[c] for c in df_cod['Codon']]
    # Calcular RSCU e relative_adaptive_weights
    aa_groups = df_cod.groupby('Amino_Acid')
    rscu = []
    rel_weights = []
    for aa in df_cod['Amino_Acid']:
        group = aa_groups.get_group(aa)
        obs = group['Obs_Freq'].values
        mean = obs.mean() if obs.mean() > 0 else 1
        rscu_val = df_cod.loc[df_cod['Amino_Acid'] == aa, 'Obs_Freq'] / mean
        max_rscu = rscu_val.max() if rscu_val.max() > 0 else 1
        rel_weight = rscu_val / max_rscu
        rscu.extend(rscu_val.tolist())
        rel_weights.extend(rel_weight.tolist())
    # Ajustar para manter a ordem dos códons válidos
    rscu_vec = [rscu[i] for i in range(len(valid_codons))]
    rel_weights_vec = [rel_weights[i] for i in range(len(valid_codons))]
    return freq_vec + rscu_vec + rel_weights_vec
    # return rscu_vec

codon_cols = [f'codon_{c}' for c in valid_codons]
rscu_cols = [f'RSCU_{c}' for c in valid_codons]
relw_cols = [f'relweight_{c}' for c in valid_codons]

######################################################

# def seq_to_kmers(seq, k=3):
# # Split sequence into k-mers.
#     seq = seq.upper().replace('N', '')  # Remove N
#     return [seq[i:i+k] for i in range(len(seq)-k+1)]

def seq_to_kmers(seq: str, k: int = 6) -> str:
    seq = seq.strip().upper().replace(" ", "")
    # gera kmers sobrepostos e junta com espaço (formato esperado pelo DNABERT v1)
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    return " ".join(kmers)

def dnabert6_embed(seqs, model, tokenizer, device, pooling="mean"):
    """
    pooling:
      - "cls": usa o embedding do token [CLS]
      - "mean": mean pooling dos tokens (ignorando padding)
    retorna: tensor [batch, hidden_size]
    """
    # if device is None:
    #     device = "cuda" if torch.cuda.is_available() else "cpu"

    # tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    # model = AutoModel.from_pretrained(model_name).to(device)
    # model.eval()

    # DNABERT-6 espera kmers separados por espaço
    kmers_text = [seq_to_kmers(s, k=6) for s in seqs]

# 1) tokeniza (pode retornar listas)
    enc = tokenizer(kmers_text, add_special_tokens=True, truncation=True)

    # 2) pad + converte para tensores torch
    enc = tokenizer.pad(enc, padding=True, return_tensors="pt")

    # 3) move para device
    enc = {k: v.to(device) for k, v in enc.items()}


    out = model(**enc)  # out.last_hidden_state: [B, T, H]
    hidden = out.last_hidden_state

    if pooling == "cls":
        emb = hidden[:, 0, :]  # token CLS
    elif pooling == "mean":
        # mean pooling com máscara para ignorar padding
        mask = enc["attention_mask"].unsqueeze(-1)  # [B, T, 1]
        summed = (hidden * mask).sum(dim=1)
        counts = mask.sum(dim=1).clamp(min=1)
        emb = summed / counts
    else:
        raise ValueError("pooling deve ser 'cls' ou 'mean'")

    return emb.cpu().detach().numpy()

def embedding_power_spectrum_mean(embedding):
    # Garantir que é 1D
    if embedding.ndim > 1:
        embedding = embedding.flatten()
    
    # Calcular FFT no vetor 1D
    fft_vals = np.fft.fft(embedding)
    
    # Calcular espectro de potência
    power_spectrum = np.abs(fft_vals) ** 2
    
    # Retornar como array 1D (achatar se necessário)
    return power_spectrum.flatten()

MODEL_ID = "zhihan1996/DNA_bert_6"

def load_dnabert6(model_id=MODEL_ID, device=None, hf_token=None):
    """
    Baixa/carrega tokenizer+model.
    Se hf_token for necessário, passe via argumento ou set HF_TOKEN no ambiente.
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    if hf_token is None:
        hf_token = os.environ.get("HF_TOKEN")

    tokenizer = AutoTokenizer.from_pretrained(
        model_id,
        do_lower_case=False,
        token=hf_token,          # ok mesmo se None
    )
    model = AutoModel.from_pretrained(
        model_id,
        token=hf_token,          # ok mesmo se None
    ).to(device)
    model.eval()
    return tokenizer, model, device

#############################
## TESTE DE FUSÃO DE GENES ##
#############################

def calcular_sobreposicao(coord1_inicio, coord1_fim, coord2_inicio, coord2_fim):
# Calculate overlap percentage between two coordinate ranges.
    """
    Calcula a porcentagem de sobreposição entre duas faixas de coordenadas.

    Args:
        coord1_inicio (int): Início da primeira coordenada.
        coord1_fim (int): Fim da primeira coordenada.
        coord2_inicio (int): Início da segunda coordenada.
        coord2_fim (int): Fim da segunda coordenada.

    Returns:
        float: Porcentagem de sobreposição (0.0 a 1.0).
    """
    # Determina o início e o fim da sobreposição
    inicio_sobreposicao = max(coord1_inicio, coord2_inicio)
    fim_sobreposicao = min(coord1_fim, coord2_fim)

    # Calcula o tamanho da sobreposição
    tamanho_sobreposicao = max(0, fim_sobreposicao - inicio_sobreposicao)

    # Calcula o tamanho total das coordenadas
    tamanho_coord1 = coord1_fim - coord1_inicio
    tamanho_coord2 = coord2_fim - coord2_inicio

    # Evita divisão por zero se uma coordenada tiver tamanho zero
    if tamanho_coord1 == 0 or tamanho_coord2 == 0:
        return 0.0

    # Calcula a porcentagem de sobreposição em relação à coordenada menor
    # ou à interseção sobre a união (Jaccard Index)
    # Aqui, vamos considerar a sobreposição em relação ao menor segmento
    # para a condição de "maior que 50% de sobreposição" ser mais robusta.
    # Se você quiser sobreposição Jaccard (interseção / união), me avise.
    
    # Opção 1: Porcentagem de sobreposição em relação ao comprimento da menor coordenada
    porcentagem_coord1 = tamanho_sobreposicao / tamanho_coord1
    porcentagem_coord2 = tamanho_sobreposicao / tamanho_coord2
    
    return max(porcentagem_coord1, porcentagem_coord2) # Retorna a maior porcentagem de sobreposição relativa


################################
### Help and argument parser ###
################################

arg_parser = argparse.ArgumentParser(description = "This script predict genes with different tools and create a consensus prediction")

arg_parser.add_argument("-i", "--input", help = "Phage genome in fasta format")
arg_parser.add_argument("--consensus", help= "Number of tools for initial consensus -- default: 3", default=3)
arg_parser.add_argument("--ident", help= "minimun identity for blast search (psiblast) -- default: 30", default=30)
arg_parser.add_argument("--cov", help= "minimun coverage for blast search (psiblast) -- default: 40", default=40)
arg_parser.add_argument("--output_path", help = "Path fot output folder -- Default: uses current directory", default=os.getcwd())
arg_parser.add_argument("--version", action="store_true", help="Show version information")
arg_parser.add_argument("--create_gff", action="store_true", help="Create GFF files from predictions")
arg_parser.add_argument("--create_blastdb", action="store_true", help="Create BLAST database")
arg_parser.add_argument("--fill", action="store_true", help="Force filling of missing genes empty coordinates")
arg_parser.add_argument("--threads", type=int, default=4, help="Number of threads to use -- default: 4")

if len(sys.argv)==1:
    arg_parser.print_help(sys.stdout)
    sys.exit(0) 

def show_version_info():
# Display version information for PELICAN.
    """Display version information from both __init__.py and setup.py"""
    print("=" * 60)
    print("PELICAN Version Information")
    print("=" * 60)
    
    # Get version from __init__.py
    try:
        from pelican import __version__ as init_version
        print(f"Version: {init_version}")
    except ImportError:
        print("Version: Unable to import")
    
    print("=" * 60)

def run_analysis():
# Main analysis function for PELICAN pipeline.
    """Main analysis function that contains all the PELICAN logic."""
    
    #######################
    ## Parsing arguments ##
    #######################

    args = arg_parser.parse_args()
    args_dict = vars(arg_parser.parse_args())

    # Check if version flag was passed
    if args_dict.get("version"):
        show_version_info()
        sys.exit(0)

    # Check if input file is provided when not using --version
    input_file = args_dict["input"]
    if not input_file:
        print("Error: Input file is required. Use -i/--input to specify the phage genome file.")
        arg_parser.print_help(sys.stderr)
        sys.exit(1)

    consensus = int(args_dict["consensus"])
    ident = int(args_dict["ident"])
    cov = int(args_dict["cov"])
    outdir = args_dict["output_path"]
    gff = args_dict["create_gff"]
    create_blastdb_flag = args_dict["create_blastdb"]
    fill = args_dict["fill"]
    threads = int(args_dict["threads"])

    torch.set_num_threads(threads)

    ##################################
    ### Setting starting variables ###
    ##################################
    
    # Display version information at start of analysis
    print("\n" + "=" * 60)
    print("PELICAN - Starting Analysis")
    try:
        from pelican import __version__ as init_version
        print(f"Version: {init_version}")
    except ImportError:
        print("Version: Unable to determine")
    print("=" * 60 + "\n")
    
    # getting primary script full path
    if '__file__' in globals():
        path = os.path.realpath(__file__)
    else:
        path = os.path.realpath(sys.argv[0])
    
    # getting primary script directory full path
    primary_script_path = path.split("/")[:-1]
    primary_script_path = "/".join(primary_script_path)
    
    #Check if input file existis
    check_file_exists(input_file)
    input_file = os.path.abspath(input_file)
    
    # Setup BLAST+ if not available
    # setup_blast()
    
    #Check if outdir exists
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    #Create folder for final outputs
    final_outdir = os.path.join(outdir, "consensus_results")
    if not os.path.exists(final_outdir):
         os.mkdir(final_outdir)
    
    #Create tmp folder
    tmp_dir = os.path.join(outdir, "tmp_dir")
    if not os.path.exists(tmp_dir):
         os.mkdir(tmp_dir)
    
    #Databases dir
    databases = os.path.join(primary_script_path, "databases")
    models_dir = os.path.join(primary_script_path, "models")
    
    # Ensure databases are uncompressed
    phrogs_db = os.path.join(databases, "Phrogs_genes.faa")
    cogs_db = os.path.join(databases, "COGorg24.faa")
    download_COG(cogs_db)
    ensure_database_uncompressed(cogs_db)
    ensure_database_uncompressed(phrogs_db)

    ## CONCATENATE ALL DATABASES FOR BLAST SEARCH ##
    combined_db = os.path.join(databases, "combined_db.faa")
    if not os.path.exists(combined_db):
        concat_fasta_files([phrogs_db, cogs_db], combined_db)
    
    #BLAST DB
    blast_db = os.path.join(databases, "combined_db.faa")
    if not os.path.exists(f"{blast_db}.pin"):
        make_blast_db(combined_db)
    elif create_blastdb_flag:
        print(color_text("Creating BLAST database as per user request..."))
        make_blast_db(combined_db)

    
    #######################
    ## Starting analysis ##
    #######################

    #General genome information (from fasta)
    for s in SeqIO.parse(input_file, "fasta"):
        genome = str(s.seq)
        genome_id = s.id
        size = len(genome)
    if gff:
        print(color_text(f"GFF flag detected, creating GFF folder in {tmp_dir}", "yellow"))
        gff_output = os.path.join(tmp_dir, "GFFs")
        if not os.path.exists(gff_output):
            os.mkdir(gff_output)


    #######################
    ## Running phanotate ##
    #######################
    print(color_text("Running Phanotate"))

    with subprocess.Popen(f"phanotate.py {input_file}",shell=True, stdout=subprocess.PIPE) as pred_process:
        phanotate_df = pd.read_csv(pred_process.stdout, sep="\t", comment="#", header=None,
                        names=["Start", "End", "Strand","Contig","SCORE"], index_col=False)

    #Getting predictions and extracting sequences
    for index, row in tqdm(phanotate_df.iterrows()):
        start = row["Start"]  # Phanotate coordinates are 1-based
        end = row["End"]
        strand = 1 if row["Strand"] == "+" else -1

        # Ensure start < end for consistent coordinate handling
        start, end, strand = ensure_start_less_than_end(start, end, strand)

        # Update the dataframe with corrected coordinates
        phanotate_df.loc[index, "Start"] = start
        phanotate_df.loc[index, "End"] = end
        phanotate_df.loc[index, "Strand"] = "+" if strand == 1 else "-"

        # Use helper function for consistent coordinate handling
        sequence = extract_sequence_from_genome(genome, start, end, strand)
        phanotate_df.loc[index, "Sequence"] = sequence

    #Adjusting table format

    phanotate_df["ORF"] = phanotate_df.apply(lambda x: f"{x['Contig']}_ORF_{x.name}", axis=1)
    phanotate_df["Strand"] = phanotate_df["Strand"].apply(lambda x: 1 if x == "+" else -1)
    phanotate_df = phanotate_df[["ORF", "Start", "End", "Strand", "Sequence"]]
    phanotate_df["Tool"] = "Phanotate"

       ## Creating fasta from Phanotate table
    tmp_phanotate_fasta = os.path.join(tmp_dir, "Phanotate_tmp.fasta")
    tmp_phanotate_GBK = os.path.join(tmp_dir, "Phanotate_GBK.gbk")
    
    with open(tmp_phanotate_fasta, "w") as file:
        for i,r in phanotate_df.iterrows():
            file.write(f">{r['ORF']}\n{r['Sequence']}\n")

    ## Creating GBK from fasta
    sequences = list(SeqIO.parse(tmp_phanotate_fasta, "fasta"))
    
    for seq in sequences:
      seq.annotations['molecule_type'] = 'DNA'
    
    SeqIO.write(sequences, tmp_phanotate_GBK, "genbank")

    ############################################################
    ## If the user wants GFF output create GFF from Phanotate ##
    ############################################################
    if gff:
        phanotate_gff = os.path.join(gff_output, f"Phanotate.gff")
        with open(phanotate_gff, "w") as gff:
            gff.write("##gff-version 3\n")
            for index, row in phanotate_df.iterrows():
                if row["Strand"] == 1:
                    strand="+"
                else:
                    strand="-"
                
                # Usar tipo de feature apropriado (tRNA ou CDS)
                if row["Tool"].startswith('tRNA'):
                    feature_type = "tRNA"
                    product = row["Tool"]
                else:
                    feature_type = "CDS"
                    product = "Phanotate predicted protein"

                gff.write(f"{genome_id}\t.\t{feature_type}\t{row['Start']}\t{row['End']}\t.\t{strand}\t.\tID={row['ORF']};product={product}\n")


    # No need to adjust coordinates - they are already correct

    ## Creating fasta from Phanotate table
    tmp_phanotate_fasta = os.path.join(tmp_dir, "Phanotate_tmp.fasta")
    tmp_phanotate_GBK = os.path.join(tmp_dir, "Phanotate_GBK.gbk")

    with open(tmp_phanotate_fasta, "w") as file:
        for i,r in phanotate_df.iterrows():
            file.write(f">{r['ORF']}\n{r['Sequence']}\n")

    ## Creating GBK from fasta
    sequences = list(SeqIO.parse(tmp_phanotate_fasta, "fasta"))

    for seq in sequences:
        seq.annotations['molecule_type'] = 'DNA'

    SeqIO.write(sequences, tmp_phanotate_GBK, "genbank")

###############
## GeneMarkS ##
###############
    print(color_text("Running GeneMarkS"))

    #########################################
    ## Looking for the GenemarkS container ##
    #########################################
    podman_container = os.path.join(primary_script_path, "genemarks", "podman_genemarks", "genemark_image.tar.gz")
    if not os.path.exists(podman_container):
        print(color_text("GeneMarkS container not found"))
        sys.exit(1)
    else:
        print(color_text("GeneMarkS container found"))
    
    genemark_install = subprocess.run(f"podman images | grep genemark", shell=True, text=True, stdout=subprocess.PIPE)
    if genemark_install.returncode != 0:
        print(color_text("GeneMarkS container not found in podman, importing...", "yellow"))
        import_cmd = f"podman load -i {podman_container}"
        subprocess.run(import_cmd, shell=True)
        print(color_text("GeneMarkS container imported successfully", "green"))
    else:
        print(color_text("GeneMarkS container already imported in podman", "green"))

    # # Caminho para a pasta do GeneMarkS
    genemarks_dir = os.path.abspath(os.path.join(primary_script_path, "genemarks"))
    # # print(genemarks_dir)
    # # Usar wrapper para rodar GeneMarkS-2
    try:
        # gms = GeneMarkS2Wrapper(genemarks_dir)
    #     gms_result = gms.predict_genes(
    #         sequence_file=input_file, 
    #         tmp_dir=tmp_dir,
    #         output_format="lst"
    #     )
    #     gms_lst = os.path.join(tmp_dir, "predictions.lst")

        ## TESTE COM GENEMARKS NO PODMAN ##

        genemarks_results = run_podman(input_file, tmp_dir)
        fasta_file_name = os.path.basename(input_file)
        gms_lst = os.path.join(tmp_dir,f"{fasta_file_name}.lst")
        print(gms_lst)
        genes = parse_predictions(output_file=gms_lst)
        # print(genes)
        # Log do conteúdo do arquivo .lst para depuração
        # print("[DEBUG] Conteúdo do arquivo .lst antes do parse:")
        # with open(gms_lst, "r") as f:
        #     print(f.read())
        # Montar DataFrame no mesmo formato dos outros preditores
        ### TESTE ###
        teste = os.path.join(tmp_dir, "GeneMarkS_tmp.csv")
        teste_df = pd.DataFrame(genes)
        teste_df.to_csv(teste, index=False)
        ### TESTE ###
        genemarks_parse = []
        genome = str(next(SeqIO.parse(input_file, "fasta")).seq)
        for i, gene in enumerate(genes):
            start = gene['start']
            end = gene['end']
            strand = 1 if gene['strand'] == "+" else -1
            start, end, strand = ensure_start_less_than_end(start, end, strand)
            sequence = extract_sequence_from_genome(genome, start, end, strand)
            orf_id = f"GeneMarkS_ORF_{i}"
            genemarks_parse.append([orf_id, start, end, strand, sequence])
        genemarks_df = pd.DataFrame(genemarks_parse, columns=["ORF", "Start", "End", "Strand", "Sequence"])
        genemarks_df["Tool"] = "GeneMarkS"
        # Criar fasta e GBK do GeneMarkS
        tmp_genemarks_fasta = os.path.join(tmp_dir, "GeneMarkS_tmp.fasta")
        tmp_genemarks_GBK = os.path.join(tmp_dir, "GeneMarkS_GBK.gbk")
        with open(tmp_genemarks_fasta, "w") as file:
            for i, r in genemarks_df.iterrows():
                file.write(f">{r['ORF']}\n{r['Sequence']}\n")
        sequences = list(SeqIO.parse(tmp_genemarks_fasta, "fasta"))
        for seq in sequences:
            seq.annotations['molecule_type'] = 'DNA'
        SeqIO.write(sequences, tmp_genemarks_GBK, "genbank")
    except Exception as e:
        print(color_text(f"Error on GeneMarkS: {e}", "red"))
        genemarks_df = pd.DataFrame(columns=["ORF", "Start", "End", "Strand", "Sequence", "Tool"])
    
    ############################################################
    ## If the user wants GFF output create GFF from GeneMarkS ##
    ############################################################
    if gff:
        genemarks_gff = os.path.join(gff_output, f"GeneMarkS.gff")
        with open(genemarks_gff, "w") as gff:
            gff.write("##gff-version 3\n")
            for index, row in genemarks_df.iterrows():
                if row["Strand"] == 1:
                    strand="+"
                else:
                    strand="-"
                
                # Usar tipo de feature apropriado (tRNA ou CDS)
                if row["Tool"].startswith('tRNA'):
                    feature_type = "tRNA"
                    product = row["Tool"]
                else:
                    feature_type = "CDS"
                    product = "GeneMarkS predicted protein"

                gff.write(f"{genome_id}\t.\t{feature_type}\t{row['Start']}\t{row['End']}\t.\t{strand}\t.\tID={row['ORF']};product={product}\n")

    ######################
    ## Running prodigal ##
    ######################
    if len(genome) < 20000:
        print(color_text("Running Prodigal in meta mode"))
        prodigal_run = pyrodigal_gv.ViralGeneFinder(meta=True)
        # prodigal_run.train(genome)
    if len(genome) >= 20000:
        print(color_text("Running Prodigal"))
        prodigal_run = pyrodigal_gv.ViralGeneFinder(meta=True)
        # prodigal_run.train(genome)
    
    #Getting predictions
    prodigal_results = []
    prodigal_genes = prodigal_run.find_genes(genome)
    for i, gene in enumerate(prodigal_genes):
        # Prodigal coordinates are 1-based, we keep them as 1-based for consistency
        # gene.sequence() already extracts the correct sequence
        # Ensure start < end for consistent coordinate handling
        start, end, strand = ensure_start_less_than_end(gene.begin, gene.end, gene.strand)
        prodigal_results.append([f"{genome_id}_ORF_{i}", start, end, strand, gene.sequence()])
    
    
    prodigal_df = pd.DataFrame.from_records(prodigal_results, columns=["ORF", "Start", "End", "Strand", "Sequence"])
    prodigal_df["Tool"] = "Prodigal"

    ############################################################
    ## If the user wants GFF output create GFF from ProdigalGV ##
    ############################################################
    if gff:
        prodigal_gff = os.path.join(gff_output, f"ProdigalGV.gff")
        with open(prodigal_gff, "w") as gff:
            gff.write("##gff-version 3\n")
            for index, row in prodigal_df.iterrows():
                if row["Strand"] == 1:
                    strand="+"
                else:
                    strand="-"
                
                # Usar tipo de feature apropriado (tRNA ou CDS)
                if row["Tool"].startswith('tRNA'):
                    feature_type = "tRNA"
                    product = row["Tool"]
                else:
                    feature_type = "CDS"
                    product = "ProdigalGV predicted protein"

                gff.write(f"{genome_id}\t.\t{feature_type}\t{row['Start']}\t{row['End']}\t.\t{strand}\t.\tID={row['ORF']};product={product}\n")

    
    # No need to adjust coordinates - Prodigal gives us correct 1-based coordinates
    # and gene.sequence() extracts the correct sequence automatically
    
    ## Creating fasta from prodigal table
    tmp_prodigal_fasta = os.path.join(tmp_dir, "Prodigal_tmp.fasta")
    tmp_prodigal_GBK = os.path.join(tmp_dir, "Prodigal_GBK.gbk")
    
    with open(tmp_prodigal_fasta, "w") as file:
        for i,r in prodigal_df.iterrows():
            file.write(f">{r['ORF']}\n{r['Sequence']}\n")

    ## Creating GBK from fasta
    sequences = list(SeqIO.parse(tmp_prodigal_fasta, "fasta"))
    
    for seq in sequences:
      seq.annotations['molecule_type'] = 'DNA'
    
    SeqIO.write(sequences, tmp_prodigal_GBK, "genbank")
    
    
    ####################
    ## Running PROKKA ##
    ####################
    print(color_text("Running Prokka Virus"))
    
    output_prokka = os.path.join(tmp_dir, "prokka_out")
    prefix = genome_id
    
    prokka_run = subprocess.run(f"prokka --kingdom Viruses --noanno --fast  --outdir {output_prokka} --prefix {prefix} --locustag {prefix} --force {input_file}",
                                shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #Getting predictions
    gff_out = os.path.join(output_prokka, f"{prefix}.gff")

    #################################
    ## Coping prokka gff to folder ##
    #################################
    if gff:
        prokka_gff = os.path.join(gff_output, f"ProkkaVirus.gff")
        shutil.copyfile(gff_out, prokka_gff)
    
    prokka_virus = pd.read_table(gff_out, comment="#", usecols=[3,4,6,8],names=["Start", "End", "Strand","ORF"], header=None).dropna()
    prokka_virus["ORF"] = prokka_virus["ORF"].apply(lambda x: x.split(";")[0]).apply(lambda x: x.replace("ID=",""))
    # GFF coordinates are 1-based, keep them as 1-based for consistency
    prokka_virus["Start"] = prokka_virus["Start"].astype(int)
    prokka_virus["End"] = prokka_virus["End"].astype(int)
    
    # Ensure start < end for consistent coordinate handling
    for index, row in prokka_virus.iterrows():
        start = int(row["Start"])
        end = int(row["End"])
        strand = 1 if row["Strand"] == "+" else -1
        
        start, end, strand = ensure_start_less_than_end(start, end, strand)
        
        prokka_virus.loc[index, "Start"] = start
        prokka_virus.loc[index, "End"] = end
        prokka_virus.loc[index, "Strand"] = "+" if strand == 1 else "-"
    
    #Concatenating predicted gene sequences
    ffn_out = os.path.join(output_prokka, f"{prefix}.ffn")
    
    ids_ = []
    sequences = []
    for s in SeqIO.parse(ffn_out, "fasta"):
        ids_.append(s.id)
        sequences.append(str(s.seq))
    map_sequences = dict(zip(ids_, sequences))
    prokka_virus["Sequence"] = prokka_virus["ORF"].map(map_sequences)
    prokka_virus["Tool"] = "ProkkaVirus"
    prokka_virus["Strand"] = prokka_virus["Strand"].apply(lambda x: 1 if x == "+" else -1)
    prokka_virus = prokka_virus[["ORF", "Start", "End", "Strand", "Sequence", "Tool"]]
    
    # No need to adjust coordinates - Prokka sequences are already correctly extracted
    
    ## Creating GBK from fasta
    tmp_prokkavirus_GBK = os.path.join(tmp_dir, "ProkkaVirus_GBK.gbk")
    sequences = list(SeqIO.parse(ffn_out, "fasta"))
    
    for seq in sequences:
      seq.annotations['molecule_type'] = 'DNA'
    
    SeqIO.write(sequences, tmp_prokkavirus_GBK, "genbank")

    ###########################
    ## Concatenating results ##
    ###########################
    # Inclui Phanotate no consenso
    complete_predictions = pd.concat([prodigal_df, prokka_virus, genemarks_df, phanotate_df], ignore_index=True)
    
    #######################################
    ## Creating comparison plot from GBK ##
    #######################################
    print(color_text("Starting plot comparison in BLAST mode"))
    
    gbk_files = (tmp_prodigal_GBK, tmp_prokkavirus_GBK, tmp_genemarks_GBK,tmp_phanotate_GBK)
    
    gbk_list = list(map(Genbank, gbk_files))
    
    #Start plot tracks
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar(ymargin=0.5)
    
    # Plot CDS features
    for gbk in gbk_list:
        track = gv.add_feature_track(gbk.name, gbk.get_seqid2size(), align_label=False)
        for seqid, features in gbk.get_seqid2features("CDS").items():
            segment = track.get_segment(seqid)
            segment.add_features(features, plotstyle="bigarrow", fc="yellow", lw=1)
    
    # Run BLAST alignment & filter by user-defined threshold
    align_coords = Blast(gbk_list, seqtype="nucleotide").run()
    align_coords = AlignCoord.filter(align_coords, length_thr=100, identity_thr=30)
    
    # Plot BLAST alignment links
    if len(align_coords) > 0:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        color, inverted_color = "green", "red"
        for ac in align_coords:
            #, vmin=min_ident
            gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity)
        gv.set_colorbar([color, inverted_color], vmin=0, vmax=100)
    
    output_compare_plot = os.path.join(final_outdir, "GBK_blast_comparison.png")
    gv.savefig(output_compare_plot, dpi=300)
    
    print(color_text(f"Plot done! --> {output_compare_plot}"))
    
    ####################################
    ## Gene comparison with UpsetPlot ##
    ####################################
    print(color_text("Starting UpsetPlot analysis"))
    
    upset_output = os.path.join(final_outdir, "UpsetPlot_shared_genes.png")
    
    #Creating unique ID based on gene coordinates
    complete_predictions["Gene_Unique_ID"] = complete_predictions.apply(lambda x: f"{x['Start']}-{x['End']}", axis=1)
    complete_predictions["COUNT"] = 1
    
    presence_abscence_matrix = complete_predictions.pivot_table("COUNT",index="Tool",columns="Gene_Unique_ID").fillna(0)
    
    #Creating the dataframe for the upsetplot
    upset_plot_df = presence_abscence_matrix == 1
    upset_plot_df = upset_plot_df.T
    upset_plot_df = upset_plot_df.value_counts()
    
    try:
        up_plot = upsetplot.UpSet(upset_plot_df, show_counts=True, intersection_plot_elements=7).plot()
        plt.title("Shared genes between tools", fontsize=20)
        plt.savefig(upset_output, dpi=300, facecolor="w")
    except Exception as e:
        print(color_text(f"Warning: UpsetPlot failed ({e}), skipping plot."))
    
    ############################
    ## First consensus filter ##
    ############################
    print(color_text("Getting first consensus"))
    
    first_consensus = presence_abscence_matrix == 1
    first_consensus = first_consensus.T
    
    first_consensus_filter = first_consensus[first_consensus.sum(axis=1) >= consensus]
    
    #Adding tool found in consensus
    tool_consensus_df = first_consensus_filter.copy()
    for index, row in first_consensus_filter.iterrows():
        columns = first_consensus_filter.columns
        tmp_list = []
        for c in columns:
            if row[c]:
                tmp_list.append(c)
            if len(tmp_list) == 3:
                con = "All"
            elif len(tmp_list) < 3:
                con = "/".join(tmp_list)
            tool_consensus_df.loc[index, "Tool consensus"] = con
    
    first_consensus_filter = tool_consensus_df.copy()
    
    print(color_text(f"First consensus found {len(first_consensus_filter)} genes", "yellow"))
    
    ###############################
    ## Checking overlapped genes ##
    ###############################

    # Nova lógica: todos os genes não consenso passam para as próximas etapas, sem remoção por sobreposição
    # print(color_text("Pulando checagem de sobreposição: todos os genes não consenso serão mantidos para análise estrutural", "yellow"))
    selected_genes = complete_predictions[complete_predictions["Gene_Unique_ID"].isin(first_consensus_filter.index)]
    check_overlap = complete_predictions[~complete_predictions["Gene_Unique_ID"].isin(first_consensus_filter.index)]
    kept_for_structural_analysis = check_overlap.copy()
    # Variáveis to_remove_overlapping, selected_list, itself_ist não são mais usadas
    print(color_text(f"{len(kept_for_structural_analysis)} genes kept for further analysis after first consensus", "yellow"))
    
    if kept_for_structural_analysis.empty:
        print(color_text("No genes kept for Blast and CNN. Skipping BLAST and CNN steps.", "yellow"))
        found_in_blast = []
        score_index = []
        cnn_only = []
    else:
        ####################################################
        ## Getting codon usage from first consensus genes ##
        ####################################################
        print(color_text("Calculating codon usage from first consensus genes", "yellow"))
        
        
        genes_for_codon = selected_genes[["Sequence"]].values.tolist()

        tmp_codon_usage = []

        for g in genes_for_codon:
            cod_u = codon_usage(g[0])
            # Se for DataFrame, pega a primeira linha como dict
            if isinstance(cod_u, pd.DataFrame):
                cod_u = cod_u.iloc[0].to_dict()
            elif hasattr(cod_u, 'shape') and len(cod_u.shape) == 2 and cod_u.shape[0] == 1:
                cod_u = cod_u[0]
            tmp_codon_usage.append(cod_u)
        codon_usage_df = pd.DataFrame(tmp_codon_usage)
        mean_codon_usage = codon_usage_df.mean()
        

        ####################################
        ## STARTING HMM PROBABILITY TESTE ##
        ####################################
        # Importa o script HMM
        print(color_text("Starting HMM analysis based on Rymuza et. al (https://doi.org/10.1093/nar/gkae685)", "yellow"))


        from pelican.models.hmm_consensus import build_signal_tracks, evaluate_genes_hmm

        # Treina o HMM apenas com os genes do consenso
        region_sets_train = []
        for idx, row in selected_genes.iterrows():
            region_sets_train.append([(int(row['Start']), int(row['End']))])

        genome_length = int(complete_predictions['End'].max())

        train_starts, train_overlaps, train_ends = build_signal_tracks(region_sets_train, genome_length)

        # Calcula matriz de transição usando apenas o consenso
        _, _, consensus_transmat = evaluate_genes_hmm(train_starts, train_overlaps, train_ends, region_sets_train, n_states=3, dynamic_transmat=True)

        # print(color_text(f"HMM transition matrix (consensus genes):\n{consensus_transmat}"))

        # Aplica o modelo nos genes não consenso
        region_sets = []
        for idx, row in kept_for_structural_analysis[kept_for_structural_analysis["Tool"] == "Phanotate"].iterrows():
            region_sets.append([(int(row['Start']), int(row['End']))])

        starts, overlaps, ends = build_signal_tracks(region_sets, genome_length)
        filtered_genes, hmm_score, transmat = evaluate_genes_hmm(starts, overlaps, ends, region_sets, n_states=3, dynamic_transmat=False, transmat=consensus_transmat)

        # print(color_text(f"HMM score (non-consensus): {hmm_score}"))
        # print(color_text(f"HMM transition matrix (non-consensus):\n{transmat}"))

              # Atualiza kept_for_structural_analysis para conter apenas genes filtrados
        phanotate_mask = kept_for_structural_analysis["Tool"] == "Phanotate"
        phanotate_filtered = kept_for_structural_analysis[phanotate_mask][
            kept_for_structural_analysis[phanotate_mask].apply(
                lambda row: any(
                    row['Start'] == start and row['End'] == end
                    for start, end, _ in filtered_genes
                ),
                axis=1
            )
        ]
        other_tools = kept_for_structural_analysis[~phanotate_mask]
        kept_for_structural_analysis = pd.concat([other_tools, phanotate_filtered], ignore_index=True)
        print(color_text(f"{len(kept_for_structural_analysis)} non-consensus genes retained after HMM filtering", "yellow"))

        ###################################
        ## Blast analysis against PHROGS ##
        ###################################
        print(color_text("Starting BLAST"))
        
        # phrogs_db = os.path.join(databases, "Phrogs_genes.faa")
        tmp_genes_fasta = os.path.join(tmp_dir, "genes_for_blast.fasta")
        
        with open(tmp_genes_fasta, "w") as file:
            for index, row in kept_for_structural_analysis.iterrows():
                seq_id = row["Gene_Unique_ID"]
                sequence = Seq(row["Sequence"]).translate(to_stop=True, stop_symbol="")
                file.write(f">{seq_id}\n{str(sequence)}\n")
        
        blast_run =  subprocess.run(f"psiblast -query {tmp_genes_fasta} -db {blast_db} -outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -max_hsps 5 -evalue 1e-5 -num_threads {threads}", shell=True, stdout=subprocess.PIPE, text=True)
        blast_out = pd.read_csv(io.StringIO(blast_run.stdout), sep="\t",names="qseqid sseqid pident qcovs evalue bitscore".split(" "), header=None)


        blast_out = blast_out[blast_out["evalue"] <= 1e-5]
        ## Score for blast ##
        blast_out["Score"] = blast_out.apply(lambda row: blast_score(row["pident"], row["qcovs"], row["evalue"]), axis=1)

        # Adiciona coluna Tool ao blast_out, mapeando qseqid para a ferramenta de predição
        geneid_to_tool = complete_predictions.set_index("Gene_Unique_ID")["Tool"].to_dict()
        blast_out["Tool_Blast"] = blast_out["qseqid"].map(geneid_to_tool)

        # blast_out = blast_out[(blast_out["pident"] >= ident) & (blast_out["qcovs"] >= cov)]

        ## Keeping variable "found in blast" empty so all proteins go to the autoencoder step
        blast_out = blast_out.sort_values(by="Score", ascending=False).drop_duplicates(subset="qseqid", keep="first")
        found_in_blast = []

        print(color_text(f"{len(blast_out)} hits scored using BLAST results"))

        ##############################
        ## AutoEncoder CNN analysis ##
        ##############################
        
        print(color_text("Starting CNN analysis"))
        
        # #READ MODEL
        # iso_forest_model  = joblib.load("/home/usuario/Projects/Setubis_BigD/ML_protein_validation/data/intermediate/ferramenta_com_ML/PAC_model.sav")
        
        for_foldseek = kept_for_structural_analysis[~kept_for_structural_analysis["Gene_Unique_ID"].isin(found_in_blast)]
        # #Fuse overllaped regions to create possible genes for analysis -- 1% overlap
        # fused_coords = fundir_anotacoes_por_coordenadas2(for_foldseek[["Gene_Unique_ID", "Start", "End"]], limite_sobreposicao=0.01)
        
        # #Add fused coordinates to dataframe
        # max_oindex = for_foldseek.index.max() + 1
        # for i, coords in enumerate(fused_coords):
        #     start, end = coords
        #     # Create a unique ID for the fused gene
        #     gene_unique_id = f"{start}-{end}"
        #     # Add the new row to the DataFrame
        #     for_foldseek.loc[max_oindex + i] = [f"FusedORF_{i}", start, end, 1,"","FUSED",gene_unique_id,1]  # Assuming strand is always 1 for simplicity
        for_foldseek.to_csv(os.path.join(final_outdir, "FOR_CNN.csv"))
        # tmp_proteins = os.path.join(tmp_dir, "proteins_for_FoldSeek.faa")
        
        #Get features from user genes
        
        coords_df = for_foldseek.copy()
        coords_df["Start"] = coords_df["Start"].astype(int)
        coords_df["End"] = coords_df["End"].astype(int)
        
        genome_record = next(SeqIO.parse(input_file, "fasta"))
        genome_seq = str(genome_record.seq)

        
        gene_data_list = []
        for index, row in coords_df.iterrows():
            phage_name = genome_id
            phage_orf = row["Gene_Unique_ID"]
            start = int(row['Start'])  # Already 1-based
            end = int(row['End'])      # Already 1-based
        
            
            # Calculate context region with 50 bases before and after (1-based coordinates)
            context_start = max(1, start)
            context_end = min(len(genome_seq), end)
            # Use helper function for consistent coordinate handling
            gene_seq = extract_sequence_from_genome(genome_seq, context_start, context_end, strand=1)
                
            # Criar um DataFrame com as coordenadas e a sequência do gene
            gene_data = {
                'GeneUniqueID': phage_orf,
                'Start': start,
                'End': end,
                'Gene_Sequence': gene_seq
            }
            gene_data_list.append(gene_data)
        
        forCNN_df = pd.DataFrame(gene_data_list)
        # print(forCNN_df.shape)
        #######################
        ## get DNABert model ##
        #######################

        tokenizer, model, device = load_dnabert6()
        #Input data
        scaler = joblib.load(os.path.join(models_dir, "minmax_scaler.pkl"))

        expected_cols = scaler.feature_names_in_ if hasattr(scaler, 'feature_names_in_') else power_input_df.columns
        adjusted_input = []
        for sequence in forCNN_df["Gene_Sequence"].values:
            embedding = dnabert6_embed([sequence], model, tokenizer, device)
            power_spec = embedding_power_spectrum_mean(embedding)
            adjusted_input.append(power_spec)
        power_input_df = pd.DataFrame(adjusted_input)
        power_input_df.columns = [f"Embedding_{i}" for i in range(1, power_input_df.shape[1] + 1)]
        # print("Done with EMbedding")
        ## ADD RSCU
        codon_feats = forCNN_df["Gene_Sequence"].apply(codon_features)
        codon_feats_df = pd.DataFrame(codon_feats.tolist(), columns=codon_cols + rscu_cols + relw_cols)
        # Concatenate the power spectrum features with the RSCU features
        complete_data = pd.concat([power_input_df.reset_index(drop=True), codon_feats_df.reset_index(drop=True)], axis=1)
        complete_data = complete_data[expected_cols]
        complete_data.to_csv(os.path.join(final_outdir, "Complete_CNN_input.csv"), index=False)
        # norm_input = np.array(power_input_df).astype(np.float32)
        # max_value = norm_input.max()
        # normalized = power_input_df.apply(min_max, axis=0)
        normalized = scaler.transform(complete_data)
        normalized_df = pd.DataFrame(normalized, columns=expected_cols)
        normalized_df.to_csv(os.path.join(final_outdir, "CNN_input.csv"), index=False)
        print(normalized.shape)
        
        #Prepare data for CNN
        # Force CPU usage for all torch operations
        
        input_dim = complete_data.shape[1]
        latent_dim = 128  # Dimensão do espaço latente
        cnn_model = CNNAutoencoder(input_dim, latent_dim=latent_dim)
        
        # Load model weights to CPU explicitly
        cnn_model.load_state_dict(torch.load(os.path.join(models_dir, "cnn_autoencoder_trained.pt"), map_location=device))
        cnn_model.to(device)
        cnn_model.eval()
        
        
        # Recuperar o erro de reconstrução para cada entrada e filtrar pelo threshold
        threshold = 0.03 #Erro medio calculando na validacao (90% dos erros)
        with torch.no_grad():
            X_valid_tensor = torch.tensor(normalized.astype(np.float32)).to(device)
            print(X_valid_tensor.shape)
            X_valid_recon = cnn_model(X_valid_tensor).cpu().numpy()
            mse_por_amostra = np.mean((X_valid_recon - normalized.astype(np.float32))**2, axis=1)
        
        CNN_SCORE = autcnn_score(mse_por_amostra)
        print("Done with CNN analysis")        


        CNN_results_df = pd.DataFrame.from_dict(dict(zip(forCNN_df["GeneUniqueID"], CNN_SCORE)), orient="index", columns=["CNN_Score"])
        # Adiciona coluna Tool ao CNN_results_df, mapeando o índice para a ferramenta de predição
        CNN_results_df["Tool_CNN"] = CNN_results_df.index.map(geneid_to_tool)

        ## Update blast score DF with AutoEncoder score ##
        final_score_data = pd.merge(blast_out, CNN_results_df, left_on="qseqid", right_index=True, how="outer").fillna(0)

        # Filtrar índices das amostras com erro acima do threshold
        # indices_abaixo = np.where(mse_por_amostra <= threshold)[0]
        # erros_abaixo = mse_por_amostra[indices_abaixo]
        
        # CNN_coords = forCNN_df.iloc[indices_abaixo]["GeneUniqueID"].tolist()
        print(color_text(f"{len(mse_por_amostra)} genes tested with the AutoEncoder model", "yellow"))
        

        ####################################################
        ## FINAL BLAST + CNN CONSENSUS BASED ON THE SCORE ##
        ####################################################
        final_score_output = os.path.join(final_outdir, "Blast_CNN_scores.csv")

        #Update CNN score based on tool used to predict the gene
        def adjust_cnn_score(row, size=len(genome)):
            tool = row["Tool_CNN"]
            score = row["CNN_Score"]
            if size >= 2000000: #Jumbo phage
                if tool == "Prodigal":
                    new_score = score * 1.1
            return new_score
        if size >= 2000000:
            for index, row in final_score_data.iterrows():
                new_score = adjust_cnn_score(row)
                final_score_data.loc[index, "CNN_Score"] = new_score

        final_score_data["Final_Score"] = (final_score_data["Score"] + final_score_data["CNN_Score"])/2
        final_score_data.to_csv(final_score_output)

        filtered_score = final_score_data[(final_score_data["Final_Score"] >= 0.67)]
        score_index = filtered_score["qseqid"].tolist()

        # cnn_only = final_score_data[(final_score_data["Score"] == 0) | (final_score_data["Final_Score"] < 0.66)]
        cnn_only = final_score_data[(final_score_data["Score"] == 0)]
        cnn_only = cnn_only[cnn_only["CNN_Score"] >= 0.994][["qseqid", "Tool_CNN"]].values.tolist()

    #################
    ## RUN ARAGORN ##
    #################
    print(color_text("Running Aragorn for tRNA prediction"))
    
    aragorn_out = os.path.join(tmp_dir, "aragorn_out.txt")
    predict_trna_aragorn(input_file, aragorn_out)
    
    # Regex para extrair coordenadas [...]  e nome do tRNA (c opcional)
    pattern = r'>.*?\s+(tRNA-\w+\([a-z]{3}\))\s+c?\[(\d+),(\d+)\]'
    
    with open(aragorn_out, "r") as file:
        content = file.read()
        
    # Encontrar todas as correspondências
    matches = re.findall(pattern, content)
    
    print(color_text(f"ARAGORN found {len(matches)} tRNAs", "yellow"))
    rna_data = []
    for match in matches:
        trna_name = match[0]
        start_coord = int(match[1])
        end_coord = int(match[2])
        
        # Ensure start < end for consistent coordinate handling (tRNAs are always forward strand for simplicity)
        start_coord, end_coord, strand = ensure_start_less_than_end(start_coord, end_coord, 1)
        
        rna_data.append([start_coord, end_coord, trna_name])
        print(color_text(f"  tRNA: {trna_name} at {start_coord}-{end_coord}", "yellow"))
    
    #############################
    ## Getting all coordinates ##
    #############################
    print(color_text("Finishing consensus and writing files"))
    
    #First consensus
    consensus_coords = first_consensus_filter.index.str.split("-").tolist()
    consensus_coords = [[int(float(j)) for j in i] for i in consensus_coords]
    
    #Adding tool info
    #add tool consensus
    consensus_tools = first_consensus_filter["Tool consensus"].values.tolist()
    adjusted_coords = []
    for i,j in zip(consensus_coords, consensus_tools):
        # print(i,j)
        t = i + [j]
        adjusted_coords.append(t)
    
    #BLAST + CNN score
    blast_CNN_consensus = []
    for i in score_index:
        coords = i.split("-")
        coords = [int(float(j)) for j in coords] + ["Blast/CNN score"]
        blast_CNN_consensus.append(coords)
    
    #CNN only
    CNN_only_consensus = []
    for id, tool in cnn_only:
        coords = id.split("-")
        coords = [int(float(j)) for j in coords] + [f"{tool}/AutoEncoder(CNN)"]
        CNN_only_consensus.append(coords)

    # #Recuperar targets
    # foldseek_consensus = []
    # for i in CNN_coords:
    #     coords = i.split("-")
    #     coords = [int(float(j)) for j in coords] + ["AutoEncoder(CNN)"]
    #     foldseek_consensus.append(coords)

    to_sort_coords = CNN_only_consensus + blast_CNN_consensus + adjusted_coords + rna_data

    #Sorting coords #NEW#
    # sorted_consensus = sorted(to_sort_coords , key=lambda k: [k[1], k[0]])
    sorted_consensus = sorted(to_sort_coords , key=lambda k: int(k[0]))
    # sorted_consensus = sorted(to_sort_coords , key=itemgetter(0))
    # sorted_consensus = to_sort_coords.sort()
    
    # Separar tRNAs dos outros genes
    trna_coords = [coord for coord in sorted_consensus if len(coord) == 3 and coord[2].startswith('tRNA')]
    gene_coords = [coord for coord in sorted_consensus if len(coord) == 3 and not coord[2].startswith('tRNA')]
    
    #Creating final dataframe
    tmp_df = pd.DataFrame.from_records(gene_coords, columns=["Start", "End", "Tool"])
    tmp_list = list(tmp_df.apply(lambda x: f"{x['Start']}-{x['End']}", axis=1))
    tmp_tools_list = tmp_df["Tool"].values.tolist()
    
    map_tools_dict = dict(zip(tmp_list,tmp_tools_list))
    
    consensus_df = complete_predictions[complete_predictions["Gene_Unique_ID"].isin(tmp_list)].drop_duplicates(subset="Gene_Unique_ID", keep="first")
    consensus_df["Tool"] = consensus_df["Gene_Unique_ID"].map(map_tools_dict)
    
    ## FINAL OVERLAP CHECK ##

    #Cheking overlap on itself
    final_overlap_remove = []
    itself_ist = consensus_df[["Start", "End", "Gene_Unique_ID"]].drop_duplicates(subset="Gene_Unique_ID", keep="first").values.tolist()
    for index, row in consensus_df.iterrows():
        for l in itself_ist:
        # print(row["Gene_Unique_ID"])
            if row["Gene_Unique_ID"] == l[2]:
                continue
            elif row["Start"] in range(l[0],l[1]) and row["End"] in range(l[0], l[1]+1):
                # print("Overlapping", row["Gene_Unique_ID"])
                final_overlap_remove.append(row["Gene_Unique_ID"])

    consensus_df = consensus_df[~consensus_df["Gene_Unique_ID"].isin(final_overlap_remove)]
    ######################################################
    ## Validação para genes totalmente sobrepostos     ##
    ## Remove o menor quando dois genes estão 90%     ##
    ## sobrepostos, mantendo apenas o com codon usage adequado   ##
    ######################################################
    print(color_text("Checking for completely overlapping genes"))
    
    genes_to_remove = []
    consensus_list = consensus_df[["Start", "End", "Gene_Unique_ID"]].values.tolist()
    
    # Check cada gene contra todos os outros
    # Nova lógica: para cada par sobreposto, manter o gene cuja codon_usage é mais próxima da média
    consensus_df_indexed = consensus_df.set_index("Gene_Unique_ID")
    # Para evitar remover ambos em pares múltiplos, manter registro dos já removidos
    genes_to_remove = set()
    for i, gene1 in enumerate(consensus_list):
        for j, gene2 in enumerate(consensus_list):
            if i >= j:
                continue  # Evita pares repetidos e auto-comparação
            id1 = gene1[2]
            id2 = gene2[2]
            if id1 in genes_to_remove or id2 in genes_to_remove:
                continue
            start1, end1 = int(gene1[0]), int(gene1[1])
            start2, end2 = int(gene2[0]), int(gene2[1])
            # Calcular sobreposição (percentual e em bp)
            overlap = calcular_sobreposicao(start1, end1, start2, end2)
            overlap_bp = max(0, min(end1, end2) - max(start1, start2))
            if overlap >= 0.8 or overlap_bp > 300:
                # Recuperar sequências
                seq1 = consensus_df_indexed.loc[id1, "Sequence"]
                seq2 = consensus_df_indexed.loc[id2, "Sequence"]
                # Calcular codon usage
                codon1 = codon_usage(seq1)
                codon2 = codon_usage(seq2)
                # Comparar com mean_codon_usage (usar apenas colunas em comum)
                mean_vec = np.ravel(mean_codon_usage.values)
                codon1_vec = np.ravel([codon1.get(k, 0) for k in mean_codon_usage.index])
                codon2_vec = np.ravel([codon2.get(k, 0) for k in mean_codon_usage.index])
                dist1 = euclidean(codon1_vec, mean_vec)
                dist2 = euclidean(codon2_vec, mean_vec)
                if dist1 < dist2:
                    genes_to_remove.add(id2)
                    print(color_text(f"  Removing overlapping gene {id2} (codon usage farther from mean) in favor of {id1}", "yellow"))
                else:
                    genes_to_remove.add(id1)
                    print(color_text(f"  Removing overlapping gene {id1} (codon usage farther from mean) in favor of {id2}", "yellow"))
    genes_to_remove = list(genes_to_remove)
    
    # Remover genes identificados como sobrepostos
    if genes_to_remove:
        consensus_df = consensus_df[~consensus_df["Gene_Unique_ID"].isin(genes_to_remove)]
        print(color_text(f"Removed {len(genes_to_remove)} overlapping genes", "yellow"))
    else:
        print(color_text("No overlapping genes found", "green"))

    # Criar DataFrame para tRNAs
    if trna_coords:
        trna_df = pd.DataFrame(trna_coords, columns=["Start", "End", "Tool"])
        trna_df["Gene_Unique_ID"] = trna_df.apply(lambda x: f"{x['Start']}-{x['End']}", axis=1)
        trna_df["ORF"] = trna_df.apply(lambda x: f"{genome_id}_tRNA_{x.name+1}", axis=1)
        trna_df["Strand"] = 1  # Assumindo strand positiva para tRNAs (pode ser ajustado se necessário)
        trna_df["Annotation"] = trna_df["Tool"]  # O nome do tRNA já vem da ferramenta
        
        # Extrair sequência dos tRNAs
        for index, row in trna_df.iterrows():
            start = int(row["Start"])  # Aragorn coordinates are 1-based
            end = int(row["End"])      # Aragorn coordinates are 1-based
            # Use helper function for consistent coordinate handling
            sequence = extract_sequence_from_genome(genome, start, end, strand=1)
            trna_df.loc[index, "Sequence"] = sequence
        
        # Check for overlap between tRNAs and consensus genes - new logic
        # New approach: 
        # - If a single tRNA overlaps with gene(s), remove the tRNA
        # - If 2+ consecutive tRNAs overlap with gene(s), remove the gene(s) instead
        print(color_text("Checking for overlap between tRNAs and consensus genes"))
        
        # First, identify which tRNAs overlap with genes
        overlapping_trnas = []
        for trna_index, trna_row in trna_df.iterrows():
            trna_start = int(trna_row["Start"])
            trna_end = int(trna_row["End"])
            
            overlapping_genes = []
            for gene_index, gene_row in consensus_df.iterrows():
                gene_start = int(gene_row["Start"])
                gene_end = int(gene_row["End"])
                
                # Calculate overlap percentage using existing function
                overlap = calcular_sobreposicao(trna_start, trna_end, gene_start, gene_end)
                
                if overlap > 0:
                    overlapping_genes.append(gene_index)
            
            if overlapping_genes:
                overlapping_trnas.append({
                    'trna_index': trna_index,
                    'trna_start': trna_start,
                    'trna_end': trna_end,
                    'trna_name': trna_row['Tool'],
                    'overlapping_genes': overlapping_genes
                })
        
        if not overlapping_trnas:
            print(color_text("No overlapping tRNAs found", "green"))
        else:
            # Group consecutive tRNAs by checking if they are adjacent in the sorted list
            trna_groups = []
            
            if len(overlapping_trnas) == 1:
                # Only one overlapping tRNA, create a single group
                trna_groups = [overlapping_trnas]
            else:
                # Multiple overlapping tRNAs, group them by proximity
                current_group = [overlapping_trnas[0]]
                
                for i in range(1, len(overlapping_trnas)):
                    prev_trna = overlapping_trnas[i-1]
                    curr_trna = overlapping_trnas[i]
                    
                    # Check if current tRNA is consecutive to previous (allowing some gap)
                    # Consider tRNAs consecutive if they are within 300bp of each other
                    gap_threshold = 300
                    if curr_trna['trna_start'] - prev_trna['trna_end'] <= gap_threshold:
                        current_group.append(curr_trna)
                    else:
                        trna_groups.append(current_group)
                        current_group = [curr_trna]
                
                trna_groups.append(current_group)
            
            # Apply the new logic: remove single tRNAs, remove genes for consecutive tRNAs
            trnas_to_remove = []
            genes_to_remove = []
            
            for group in trna_groups:
                if len(group) == 1:
                    # Single tRNA overlapping with gene(s) - remove the tRNA
                    trna = group[0]
                    trnas_to_remove.append(trna['trna_index'])
                    print(color_text(f"  Removing single tRNA {trna['trna_name']} at {trna['trna_start']}-{trna['trna_end']}", "yellow"))
                else:
                    # Multiple consecutive tRNAs overlapping with gene(s) - remove the gene(s)
                    all_overlapping_genes = set()
                    trna_names = []
                    for trna in group:
                        all_overlapping_genes.update(trna['overlapping_genes'])
                        trna_names.append(trna['trna_name'])
                    
                    genes_to_remove.extend(list(all_overlapping_genes))
                    print(color_text(f"  Found {len(group)} consecutive tRNAs ({', '.join(trna_names)}), removing {len(all_overlapping_genes)} overlapping gene(s)", "yellow"))
            
            # Remove tRNAs and genes
            if trnas_to_remove:
                trna_df = trna_df.drop(trnas_to_remove).reset_index(drop=True)
                print(color_text(f"Removed {len(trnas_to_remove)} single overlapping tRNAs", "yellow"))
            
            if genes_to_remove:
                consensus_df = consensus_df.drop(genes_to_remove).reset_index(drop=True)
                print(color_text(f"Removed {len(genes_to_remove)} genes that overlap with consecutive tRNAs", "yellow"))
    
    # Combinar genes e tRNAs (after applying overlap resolution logic)
    if trna_coords and not trna_df.empty:
        consensus_df = pd.concat([consensus_df, trna_df], ignore_index=True)
    

    consensus_df = consensus_df.sort_values(by="Start", ascending=True)
    consensus_df.reset_index(drop=True, inplace=True)

    ###############################################
    ## Preencher buracos do consenso com genes  ##
    ## de preditores que melhor se encaixam     ##
    ###############################################
    if fill:
        print(color_text("Looking for empty regions on the consensus", "yellow"))
        # Definir mínimo de tamanho de gap para considerar (pode ser ajustado)
        min_gap_size = 30
        # Obter lista de coordenadas ordenadas
        coords = consensus_df.sort_values(by="Start")[["Start", "End"]].values.tolist()
        gaps = []
        for i in range(len(coords) - 1):
            end_current = coords[i][1]
            start_next = coords[i+1][0]
            if start_next - end_current > min_gap_size:
                gaps.append((end_current + 1, start_next - 1))

        # Buscar candidatos em complete_predictions para cada gap
        genes_adicionados = []
        for gap_start, gap_end in gaps:
            # Candidatos que estão totalmente dentro do gap e não estão no consenso
            candidatos = complete_predictions[
                (complete_predictions["Start"] >= gap_start) &
                (complete_predictions["End"] <= gap_end) &
                (~complete_predictions["Gene_Unique_ID"].isin(consensus_df["Gene_Unique_ID"]))
            ].copy()
            if not candidatos.empty:
                # Escolher o maior gene (maior cobertura do gap)
                candidatos["coverage"] = candidatos["End"] - candidatos["Start"] + 1
                melhor = candidatos.sort_values(by="coverage", ascending=False).iloc[0]
                # Marcar origem
                melhor["Tool"] = f"GapFiller_{melhor['Tool']}"
                genes_adicionados.append(melhor)
                print(color_text(f"  Gap {gap_start}-{gap_end} filled with gene {melhor['Gene_Unique_ID']} ({melhor['Tool']})", "yellow"))

        # Adicionar genes ao consenso
        if genes_adicionados:
            consensus_df = pd.concat([consensus_df, pd.DataFrame(genes_adicionados)], ignore_index=True)
            consensus_df = consensus_df.sort_values(by="Start", ascending=True).reset_index(drop=True)
            print(color_text(f"{len(genes_adicionados)} gaps filled with genes from predictors.", "yellow"))
        else:
            print(color_text("No gaps filled (no candidate genes found in gaps).", "yellow"))

    #########################################
    ## Check gene length divisibility by 3 ##
    #########################################
    # print(color_text("Checking gene length divisibility by 3"))
    
    genes_adjusted = 0
    for index, row in consensus_df.iterrows():
        # Skip tRNAs from divisibility check (they don't need to be divisible by 3)
        if row["Tool"].startswith('tRNA'):
            continue
            
        gene_length = row["End"] - row["Start"] + 1  # +1 because coordinates are inclusive
        if gene_length % 3 != 0:
            # Calculate how many nucleotides we need to remove to make it divisible by 3
            remainder = gene_length % 3
            # Adjust start position to make length divisible by 3
            adjustment_needed = remainder
            
            # Move start position forward to reduce gene length
            old_start = row["Start"]
            new_start = old_start + adjustment_needed
            consensus_df.loc[index, "Start"] = new_start
            genes_adjusted += 1
            
            # Verify the new length is divisible by 3
            new_length = row["End"] - new_start + 1
            # print(color_text(f"Gene {row['Gene_Unique_ID']}: adjusted Start from {old_start} to {new_start} (length: {gene_length} -> {new_length}, remainder: {remainder} -> {new_length % 3})", "yellow"))
    
    if genes_adjusted > 0:
        # print(color_text(f"Total genes adjusted for divisibility by 3: {genes_adjusted}", "yellow"))
        
        # Atualizar as sequências baseadas nas novas coordenadas
        # print(color_text("Updating gene sequences based on adjusted coordinates"))
        for index, row in consensus_df.iterrows():
            # Skip tRNAs from sequence update (they already have correct sequences)
            if row["Tool"].startswith('tRNA'):
                continue
                
            start = int(row["Start"])  # 1-based coordinates
            end = int(row["End"])      # 1-based coordinates
            strand = int(row["Strand"])
            
            # Use helper function for consistent coordinate handling
            new_sequence = extract_sequence_from_genome(genome, start, end, strand)
            consensus_df.loc[index, "Sequence"] = new_sequence
    
    # Final verification that all genes are divisible by 3
    # print(color_text("Final verification of gene lengths"))
    problematic_genes = 0
    for index, row in consensus_df.iterrows():
        # Skip tRNAs
        if row["Tool"].startswith('tRNA'):
            continue
            
        gene_length = row["End"] - row["Start"] + 1  # +1 because coordinates are inclusive
        if gene_length % 3 != 0:
            problematic_genes += 1
            # print(color_text(f"WARNING: Gene {row['Gene_Unique_ID']} still not divisible by 3! Length: {gene_length} (remainder: {gene_length % 3})", "red"))
    
    #Adjust ORF names (skip tRNAs as they already have proper names)
    for index, row in consensus_df.iterrows():
        if not row["Tool"].startswith('tRNA'):
            consensus_df.loc[index,"ORF"] = f"{genome_id}_{index+1}"
    
    final_table_output = os.path.join(final_outdir, "Consensus_coords_table.tsv")
    
    ################################################
    ## Create temporary faa file for BLAST search ##
    ################################################
    
    print(color_text("Creating TMP fasta files with proteins for annotation using BLAST search", "yellow"))
    TMP_faa_out = os.path.join(tmp_dir, "TMP.faa")
    
    #Protein fasta (excluding tRNAs)
    with open(TMP_faa_out, "w") as faa:
        for number, data in enumerate(consensus_df.iterrows()):
            # Skip tRNAs for protein annotation
            if data[1]["Tool"].startswith('tRNA'):
                continue
            gene_id = f"{genome_id}_{number+1}"
            sequence = Seq(data[1]["Sequence"]).translate(stop_symbol="")
            faa.write(f">{gene_id}\n{sequence}\n")
    
    ###################################
    ## Blast analysis against PHROGS ##
    ###################################
    print(color_text("Starting final functional annotation"))
    
    # phrogs_db = os.path.join(databases, "Phrogs_genes.faa")
    phrogs_annot = os.path.join(databases, "phrog_annot_v4.tsv")
    
    #Load annotation table for mapping
    annot_table = pd.read_table(phrogs_annot, usecols=["phrog", "annot"])
    annot_table["phrog"] = annot_table["phrog"].apply(lambda x: f"phrog_{x}")
    annot_table["annot"] = annot_table["annot"].combine_first(annot_table["phrog"])
    
    annot_dict = dict(zip(annot_table["phrog"].values, annot_table["annot"].values))
    
    # tmp_genes_fasta = os.path.join(tmp_dir, "genes_for_blast.fasta")
    
    blast_run =  subprocess.run(f"blastp -query {TMP_faa_out} -subject {phrogs_db} -outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -evalue 1e-5", shell=True, stdout=subprocess.PIPE, text=True)
    blast_out = pd.read_csv(io.StringIO(blast_run.stdout), sep="\t",names="qseqid sseqid pident qcovs evalue bitscore".split(" "), header=None)
    
    blast_out = blast_out[(blast_out["pident"] >= ident) & (blast_out["qcovs"] >= cov)]
    
    blast_out.sort_values(by="pident", ascending=False).drop_duplicates(subset="qseqid", keep="first", inplace=True)
    annotated_map = dict(zip(blast_out["qseqid"].values, blast_out["sseqid"].values))
    
    consensus_df["Annotation"] = consensus_df["ORF"].map(annotated_map)
    
    #################################
    ## Write final consensus table ##
    #################################
    
    # Handle annotations: tRNAs get their Tool value, others get "Unnanotated protein" if no annotation
    for index, row in consensus_df.iterrows():
        if pd.isna(row["Annotation"]):
            if row["Tool"].startswith('tRNA'):
                consensus_df.loc[index, "Annotation"] = row["Tool"]
                consensus_df.loc[index, "Tool"] = "tRNA"
            else:
                consensus_df.loc[index, "Annotation"] = "Unnanotated protein"
    
    consensus_df.replace(annot_dict, inplace=True)
    
    #Adjusting ORF names for final consensus
    consensus_df.reset_index(drop=True, inplace=True)
    
    #Adjust ORF names (skip tRNAs as they already have proper names)
    for index, row in consensus_df.iterrows():
        if not row["Tool"].startswith('tRNA'):
            consensus_df.loc[index,"ORF"] = f"{genome_id}_{index+1}"
    
    consensus_df.to_csv(final_table_output, sep="\t", index=False)
    
    ######################################
    ## Creating final fna and faa files ##
    ######################################
    print(color_text("Creating fasta files with genes and proteins"))
    final_fna_out = os.path.join(final_outdir, "consensus_genes.fna")
    final_faa_out = os.path.join(final_outdir, "consensus_proteins.faa")
    
    #Nucleotide fasta (including both genes and tRNAs)
    with open(final_fna_out, "w") as fna:
        for number, data in enumerate(consensus_df.iterrows()):
            gene_id = data[1]["ORF"]  # Use the proper ORF name (includes tRNA names)
            sequence = data[1]["Sequence"]
            fna.write(f">{gene_id}\n{sequence}\n")
    
    #Protein fasta (excluding tRNAs)
    with open(final_faa_out, "w") as faa:
        for number, data in enumerate(consensus_df.iterrows()):
            # Skip tRNAs for protein file
            if data[1]["Tool"].startswith('tRNA'):
                continue
            gene_id = data[1]["ORF"]  # Use the proper ORF name
            sequence = Seq(data[1]["Sequence"]).translate(stop_symbol="")
            faa.write(f">{gene_id}\n{sequence}\n")
    
    
    #################################
    ## Creating consensus GFF file ##
    #################################
    
    output_gff = os.path.join(final_outdir, f"{genome_id}_consensus.gff")
    
    with open(output_gff, "w") as gff:
        gff.write("##gff-version 3\n")
        for index, row in consensus_df.iterrows():
            if row["Strand"] == 1:
                strand="+"
            else:
                strand="-"
            
            # Usar tipo de feature apropriado (tRNA ou CDS)
            if row["Tool"].startswith('tRNA'):
                feature_type = "tRNA"
                product = row["Tool"]
            else:
                feature_type = "CDS"
                product = row['Annotation']

            gff.write(f"{genome_id}\t.\t{feature_type}\t{row['Start']}\t{row['End']}\t.\t{strand}\t.\tID={row['ORF']};product={product}\n")

    ##############################
    ## From the GFF crete a GBK ##
    ##############################
    
    convert_gff_to_gbk(output_gff, input_file, consensus_df)
    #https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py
    
    ##################################
    ## Plot final consenus sequence ##
    ##################################
    
    print(color_text("Plotting final consensus"))
    
    output_final_consensus_plot = os.path.join(final_outdir, "Consensus_plot.png")
    output_final_consensus_plotSVG = os.path.join(final_outdir, "Consensus_plot.svg")
    
    #Get all unique tools to create random colors
    unique_tools = consensus_df["Tool"].unique().tolist() 
    
    number_of_colors = len(unique_tools)
    
    # color = mcp.gen_color(cmap="tab20",n=10)
    color = ['#3182bd', '#e6550d', '#fd8d3c', '#31a354', '#756bb1', '#9e9ac8', '#636363', '#969696', '#bdbdbd', '#d9d9d9']
    color_ = color[:number_of_colors]
    
    color_dict = dict(zip(unique_tools,color_))
    color_dict["All"] = "#12120f" #Black
    color_dict["tRNA"] = "#ff0000" #Red
    
    # plot_consensus
    gv = GenomeViz()
    gv.set_scale_xticks(ymargin=0.5)
    
    track = gv.add_feature_track(f"{genome_id} Consensus", size)
    track.add_sublabel()
    
    for i,r in consensus_df.iterrows():
        track.add_feature(r["Start"], r["End"], r["Strand"], fc=color_dict[r["Tool"]])
    
    # Create the plot
    fig = gv.plotfig()
    # Add legend using matplotlib
    legend_patches = [mpatches.Patch(color=color, label=tool) for tool, color in color_dict.items()]
    plt.legend(handles=legend_patches, loc='upper right', bbox_to_anchor=(1.15, 1))

    plt.tight_layout()
    plt.savefig(output_final_consensus_plot, dpi=300)
    plt.savefig(output_final_consensus_plotSVG)
    
    print(color_text("ALL DONE!"))


def run_pelican(input_file, consensus=2, ident=30, cov=60, output_path=None):
# Run PELICAN programmatically with given parameters.
    """
    Run PELICAN analysis programmatically.
    
    Args:
        input_file (str): Path to phage genome in fasta format
        consensus (int): Number of tools for initial consensus (default: 2)
        ident (int): Minimum identity for blast search (default: 50)
        cov (int): Minimum coverage for blast search (default: 60)
        output_path (str): Path for output folder (default: current directory)
    """
    # Override sys.argv to simulate command line arguments
    original_argv = sys.argv
    sys.argv = [
        'pelican',
        '--input', input_file,
        '--consensus', str(consensus),
        '--ident', str(ident),
        '--cov', str(cov),
        '--output_path', output_path if output_path else os.getcwd()
    ]
    
    try:
        run_analysis()
    finally:
        # Restore original sys.argv
        sys.argv = original_argv


def main():
# Main entry point for command line interface.
    """Main entry point for command line interface."""
    run_analysis()


if __name__ == "__main__":
    main()