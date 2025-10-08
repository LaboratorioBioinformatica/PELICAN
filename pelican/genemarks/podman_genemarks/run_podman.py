import subprocess
import os

# input_path = "/home/usuario/Projects/GeneMarkS/teste/phage_Aloeri.fasta"
# output_dir = "/home/usuario/Projects/GeneMarkS/teste/teste_2"

def run_podman(input_path, output_dir):
    if not os.path.isfile(input_path):
        print(f"Arquivo FASTA não encontrado: {input_path}")
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
        "-v", f"{os.path.abspath(output_dir)}:/output",  # usa caminho absoluto
        "-v", f"{key}:/gm_key.gz",
        "-w", "/opt/gmsuite",
        "genemark",
        "bash", "-c",
        (
            f"cp /data/{fasta} . && "
            f"gunzip -c /gm_key.gz > ~/.gm_key && "
            f"./gmsn.pl --phage {fasta} && "
            f"cp {fasta}.lst /output/"
        )
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # if result.returncode != 0:
    #     print("Erro ao rodar GeneMark:", result.stderr)
    # else:
    #     if os.path.exists(output_lst):
    #         with open(output_lst) as f:
    #             output = f.read()
    #         print("Resultado do GeneMark:")
    #         print(output)
    #     else:
    #         print("Arquivo de saída não encontrado:", output_lst)

# run_podman(input_path, output_dir)