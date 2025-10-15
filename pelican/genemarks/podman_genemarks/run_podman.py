import subprocess
import os

def run_podman(input_path, output_dir):
    if not os.path.isfile(input_path):
        print(f"FASTA file not found: {input_path}")
        return

    os.makedirs(output_dir, exist_ok=True)  # Garante que o diretório de saída existe

    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)

    fasta_dir = os.path.dirname(input_path)
    key = os.path.join(current_dir, "gm_key.gz")
    fasta = os.path.basename(input_path)
    output_lst = os.path.join(output_dir, fasta + ".lst")

    # Remove old files
    if os.path.exists(output_lst):
        os.remove(output_lst)

    cmd = [
        "podman", "run", "--rm", "--entrypoint", "",
        "-v", f"{os.path.abspath(fasta_dir)}:/data",
        "-v", f"{os.path.abspath(output_dir)}:/output",
        "-v", f"{key}:/gm_key.gz",
        "-w", "/opt/gmsuite",
        "genemark",
        "bash", "-c",
        (
            f"exec > /output/genemark_execution.log 2>&1 && "
            f"export LC_ALL=C && "  # Força idioma inglês
            f"echo 'Starting GeneMark processing...' && "
            f"echo 'Input file: {fasta}' && "
            f"echo 'Working directory: /opt/gmsuite' && "
            f"cp /data/{fasta} . && "
            f"echo 'FASTA file copied successfully' && "
            f"gunzip -c /gm_key.gz > ~/.gm_key && "
            f"echo 'License key extracted successfully' && "
            f"echo 'Running GeneMark...' && "
            f"./gmsn.pl --phage {fasta} && "
            f"echo 'GeneMark execution completed' && "
            f"cp {fasta}.lst /output/ && "
            f"echo 'Output file copied to /output/' && "
            f"echo 'Process finished successfully'"
        )
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)