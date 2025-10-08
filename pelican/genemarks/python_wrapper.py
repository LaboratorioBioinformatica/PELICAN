#!/usr/bin/env python3
"""
Wrapper Python para GeneMarkS-2
ATENÇÃO: Respeite os termos da licença do GeneMarkS-2
"""

import subprocess
import os
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any
import shutil

class GeneMarkS2Wrapper:
    """
    Wrapper Python para o preditor de genes GeneMarkS-2
    
    IMPORTANTE: Este wrapper requer licença válida do GeneMarkS-2
    e só pode ser usado conforme os termos da licença original.
    """
    
    def __init__(self, gms2_path: str):
        """
        Inicializa o wrapper
        
        Args:
            gms2_path: Caminho para o diretório do GeneMarkS-2
        """
        self.gms2_path = Path(gms2_path)
        self.gms2_script = self.gms2_path / "gms2.pl"
        
        # Verificar se os arquivos necessários existem
        if not self.gms2_script.exists():
            raise FileNotFoundError(f"Script gms2.pl não encontrado em {self.gms2_script}")
        
        # Verificar se a chave está instalada
        key_file = Path.home() / ".gmhmmp2_key"
        if not key_file.exists():
            raise FileNotFoundError(
                "Chave gmhmmp2_key não encontrada. "
                "Execute: gunzip -c gm_key.gz > ~/.gmhmmp2_key"
            )

    def predict_genes(self, tmp_dir: str,
                     sequence_file: str,
                     genome_type: str = "auto",
                     genetic_code: str = "auto",
                     output_file: Optional[str] = None,
                     output_format: str = "lst",
                     get_fasta: bool = False) -> Dict[str, Any]:
        """
        Executa predição de genes com GeneMarkS-2
        
        Args:
            sequence_file: Arquivo FASTA com sequência genômica
            genome_type: Tipo do genoma - opções válidas: "bacteria", "archaea", "auto" (padrão: "auto")
            genetic_code: Código genético (auto, 11, 4, 25, 15)
            output_file: Arquivo de saída (opcional)
            output_format: Formato de saída (lst, gff, gtf, gff3)
            get_fasta: Se True, também gera arquivos FASTA das predições
            
        Returns:
            Dict com resultados e caminhos dos arquivos gerados
        """
        
        # Validar genome_type
        valid_genome_types = ["bacteria", "archaea", "auto"]
        if genome_type not in valid_genome_types:
            raise ValueError(f"genome_type deve ser um de: {valid_genome_types}")
        
        # Copiar o arquivo de entrada para o diretório do GeneMarkS
        import shutil
        seq_basename = os.path.basename(sequence_file)
        seq_target = str(self.gms2_path / seq_basename)
        shutil.copy(sequence_file, seq_target)
        # # DEBUG: Listar arquivos no diretório do GeneMarkS
        # print("[GeneMarkS2Wrapper] Arquivos no diretório do GeneMarkS:")
        # for f in os.listdir(self.gms2_path):
        #     print("  ", f)
        # # DEBUG: Mostrar conteúdo do arquivo FASTA copiado
        # print(f"[GeneMarkS2Wrapper] Conteúdo do arquivo FASTA copiado ({seq_target}):")
        # with open(seq_target, 'r') as f:
        #     print(f.read())
        # Criar arquivo de saída temporário se não especificado
        if output_file is None:
            output_file = os.path.join(tmp_dir, "predictions.lst")
        output_basename = os.path.basename(output_file)
        output_target = str(self.gms2_path / output_basename)
        # Construir comando
        cmd = [
            "perl", str(self.gms2_script),
            "--seq", seq_basename,
            "--genome-type", genome_type,
            "--gcode", genetic_code,
            "--output", output_basename,
            "--format", output_format
        ]
        
        # Adicionar opções para arquivos FASTA
        fnn_file = None
        faa_file = None
        if get_fasta:
            fnn_file = output_basename.replace('.lst', '_genes.fnn')
            faa_file = output_basename.replace('.lst', '_proteins.faa')
            cmd.extend(["--fnn", fnn_file, "--faa", faa_file])
        
        try:
            # Garantir que o PATH inclua o diretório do GeneMarkS
            env = os.environ.copy()
            env["PATH"] = str(self.gms2_path) + os.pathsep + env.get("PATH", "")
            # Executar GeneMarkS-2
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.gms2_path,
                timeout=3600,  # Timeout de 1 hora
                env=env
            )
            # DEBUG: Salvar stderr do GeneMarkS-2 em arquivo de log
            log_path = os.path.join(self.gms2_path, 'genemarks_stderr.log')
            with open(log_path, 'a') as logf:
                logf.write(f"\n[GeneMarkS2Wrapper] Comando: {' '.join(cmd)}\n")
                logf.write(result.stderr)
            # Verificar se a execução foi bem-sucedida
            if result.returncode != 0:
                raise RuntimeError(f"GeneMarkS-2 falhou: {result.stderr}")
            # Ler resultados
            predictions = []
            # Copiar o arquivo de saída de volta, se gerado, e apagar da pasta do genemarks
            if os.path.exists(output_target):
                shutil.copy(output_target, output_file)
                with open(output_file, 'r') as f:
                    predictions = f.readlines()
                os.remove(output_target)
            results = {
                "success": True,
                "output_file": output_file,
                "predictions": predictions,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": " ".join(cmd)
            }
            if get_fasta:
                fnn_target = str(self.gms2_path / fnn_file) if fnn_file else None
                faa_target = str(self.gms2_path / faa_file) if faa_file else None
                if fnn_target and os.path.exists(fnn_target):
                    dest_fnn = os.path.join(os.path.dirname(output_file), os.path.basename(fnn_target))
                    shutil.copy(fnn_target, dest_fnn)
                    results["fnn_file"] = dest_fnn
                    os.remove(fnn_target)
                if faa_target and os.path.exists(faa_target):
                    dest_faa = os.path.join(os.path.dirname(output_file), os.path.basename(faa_target))
                    shutil.copy(faa_target, dest_faa)
                    results["faa_file"] = dest_faa
                    os.remove(faa_target)
            # Limpar arquivo de entrada temporário
            if os.path.exists(seq_target):
                os.remove(seq_target)
            return results
        except subprocess.TimeoutExpired:
            raise RuntimeError("GeneMarkS-2 excedeu tempo limite de execução")
        except Exception as e:
            raise RuntimeError(f"Erro ao executar GeneMarkS-2: {str(e)}")
    
    def parse_predictions(self, output_file: str) -> list:
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

# Exemplo de uso
def main():
    """Exemplo de uso do wrapper"""
    
    # Inicializar wrapper
    gms2_path = "/home/usuario/Projects/GeneMarkS2/gms2_linux_64"
    wrapper = GeneMarkS2Wrapper(gms2_path)
    
    # Criar arquivo de teste
    test_fasta = "/tmp/test_genome.fasta"
    with open(test_fasta, 'w') as f:
        f.write(">test_genome\n")
        f.write("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA\n")
        f.write("ATGCGCAAATTAAATAAAAAACACCCTTTTATGATCTGCCAACTTTAAATCGGTGGGATACGGTAC\n")
    
    try:
        # Executar predição
        results = wrapper.predict_genes(
            sequence_file=test_fasta,
            genome_type="bacteria",
            get_fasta=True
        )
        
        print(f"Predição executada com sucesso!")
        print(f"Arquivo de saída: {results['output_file']}")
        print(f"Número de linhas de predição: {len(results['predictions'])}")
        
        # Fazer parsing dos resultados
        if os.path.exists(results['output_file']):
            genes = wrapper.parse_predictions(results['output_file'])
            print(f"Genes preditos: {len(genes)}")
            for gene in genes:
                print(f"  Gene {gene['gene_id']}: {gene['start']}-{gene['end']} ({gene['strand']})")
    
    except Exception as e:
        print(f"Erro: {e}")
    
    finally:
        # Limpar arquivo de teste
        if os.path.exists(test_fasta):
            os.remove(test_fasta)

if __name__ == "__main__":
    main()
