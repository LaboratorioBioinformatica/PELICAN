"""
GeneMarkS runner for PELICAN.

Runs `gmsn.pl --phage` (GeneMarkS 4.32) with one of three backends:

- native : the same binaries shipped in the container image, copied to
           pelican/genemarks/gmsuite and executed with the local perl.
           No podman/docker/sudo needed (Linux x86-64 only).
- podman : the original container image (genemark_image.tar.gz).
- docker : the same container image, run with docker.

Backend selection: environment variable PELICAN_GENEMARK_BACKEND
(auto | native | podman | docker). Default "auto" tries native, then podman,
then docker.

License key lookup (GeneMark keys expire after 400 days):
  1. PELICAN_GM_KEY environment variable (path to gm_key, plain or .gz)
  2. ~/.gm_key (standard GeneMark location)
  3. key bundled with PELICAN (gm_key.gz in this folder)

Outputs are the same for every backend:
  <output_dir>/<fasta>.lst and <output_dir>/genemark_execution.log
"""
import gzip
import os
import platform
import shutil
import subprocess
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
GMSUITE_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "gmsuite"))
SELFTEST_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "selftest"))
IMAGE_ARCHIVE = os.path.join(CURRENT_DIR, "genemark_image.tar.gz")
BUNDLED_KEY = os.path.join(CURRENT_DIR, "gm_key.gz")
IMAGE_NAME = "genemark"
BACKENDS = ("native", "podman", "docker")
LOG_NAME = "genemark_execution.log"

# Binaries/scripts executed by gmsn.pl in --phage mode
NATIVE_EXECUTABLES = ("gmsn.pl", "gm", "gmhmmp", "probuild", "Gibbs3", "gibbs", "gc", "mkmat", "viewmat", "matinfo")

_selected_backend = None


def _color(text, color="red"):
    if not sys.stdout.isatty():
        return text
    codes = {"red": "\033[91m", "green": "\033[92m", "yellow": "\033[93m"}
    return f"{codes.get(color, '')}{text}\033[0m"


def find_gm_key():
    """Return the path of the GeneMark license key to use."""
    env_key = os.environ.get("PELICAN_GM_KEY")
    if env_key:
        if os.path.isfile(env_key):
            return os.path.abspath(env_key)
        print(_color(f"PELICAN_GM_KEY points to a missing file: {env_key}", "yellow"))
    home_key = os.path.expanduser("~/.gm_key")
    if os.path.isfile(home_key):
        return home_key
    return BUNDLED_KEY


def _read_key(key_path):
    """Read a key file, decompressing it if it is gzipped."""
    with open(key_path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(key_path, "rb") as fh:
            return fh.read()
    with open(key_path, "rb") as fh:
        return fh.read()


def license_expired(log_path):
    """True if the GeneMark log reports an expired license."""
    try:
        with open(log_path, "r", errors="ignore") as fh:
            return "license period has ended" in fh.read()
    except OSError:
        return False


def explain_license_error():
    print(_color(
        "GeneMark license key has expired (GeneMark keys are valid for 400 days).\n"
        "  1. Get a free academic key at https://exon.gatech.edu/GeneMark/license_download.cgi\n"
        "     (choose GeneMarkS / LINUX 64 and download the key file gm_key_64.gz)\n"
        "  2. Install it with:  pelican-setup --gm-key /path/to/gm_key_64.gz\n"
        "     (or: gunzip -c gm_key_64.gz > ~/.gm_key)", "red"))


###################
## Native backend ##
###################

def _ensure_executable():
    for name in NATIVE_EXECUTABLES:
        path = os.path.join(GMSUITE_DIR, name)
        if os.path.isfile(path) and not os.access(path, os.X_OK):
            try:
                os.chmod(path, os.stat(path).st_mode | 0o755)
            except OSError:
                pass


def native_available():
    """Check if the native backend can run on this machine."""
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        return False, "native GeneMarkS binaries require Linux x86-64"
    if not os.path.isfile(os.path.join(GMSUITE_DIR, "gmsn.pl")):
        return False, f"GeneMarkS files not found in {GMSUITE_DIR}"
    if shutil.which("perl") is None:
        return False, "perl not found in PATH"
    _ensure_executable()
    try:
        # gmhmmp without arguments prints its usage; a loader error means missing libraries
        probe = subprocess.run([os.path.join(GMSUITE_DIR, "gmhmmp")], capture_output=True, text=True, timeout=30)
        if "error while loading shared libraries" in (probe.stderr or ""):
            return False, probe.stderr.strip()
    except OSError as e:
        return False, f"cannot execute GeneMarkS binaries: {e}"
    return True, "ok"


def _run_native(input_path, output_dir):
    fasta = os.path.basename(input_path)
    log_path = os.path.join(output_dir, LOG_NAME)
    with tempfile.TemporaryDirectory(prefix="pelican_gms_") as tmp:
        work_dir = os.path.join(tmp, "work")
        home_dir = os.path.join(tmp, "home")
        os.makedirs(work_dir)
        os.makedirs(home_dir)
        # Private HOME holding the license key: the user's home is not touched
        with open(os.path.join(home_dir, ".gm_key"), "wb") as fh:
            fh.write(_read_key(find_gm_key()))
        shutil.copy(input_path, os.path.join(work_dir, fasta))

        env = os.environ.copy()
        env["HOME"] = home_dir
        env["LC_ALL"] = "C"
        cmd = ["perl", os.path.join(GMSUITE_DIR, "gmsn.pl"), "--phage", fasta]
        with open(log_path, "w") as log:
            log.write("Starting GeneMark processing (native backend)...\n")
            log.write(f"Input file: {fasta}\n")
            log.flush()
            result = subprocess.run(cmd, cwd=work_dir, env=env, stdout=log, stderr=subprocess.STDOUT)
            produced = os.path.join(work_dir, f"{fasta}.lst")
            if result.returncode == 0 and os.path.exists(produced):
                shutil.copy(produced, os.path.join(output_dir, f"{fasta}.lst"))
                log.write("Process finished successfully\n")
    return result.returncode


######################
## Container backends ##
######################

def _selinux_enforcing():
    try:
        with open("/sys/fs/selinux/enforce") as fh:
            return fh.read().strip() == "1"
    except OSError:
        return False


def container_available(engine):
    if shutil.which(engine) is None:
        return False, f"{engine} not found in PATH"
    probe = subprocess.run([engine, "info"], capture_output=True, text=True)
    if probe.returncode != 0:
        return False, f"'{engine} info' failed: {probe.stderr.strip()[:300]}"
    return True, "ok"


def ensure_image(engine):
    """Load the GeneMarkS image into podman/docker if it is not there yet."""
    if engine == "podman":
        exists = subprocess.run([engine, "image", "exists", IMAGE_NAME]).returncode == 0
    else:
        exists = subprocess.run([engine, "image", "inspect", f"localhost/{IMAGE_NAME}:latest"],
                                capture_output=True).returncode == 0
    if exists:
        print(_color(f"GeneMarkS container already imported in {engine}", "green"))
        return True
    if not os.path.exists(IMAGE_ARCHIVE):
        print(_color(f"GeneMarkS container image not found: {IMAGE_ARCHIVE}"))
        return False
    print(_color(f"GeneMarkS container not found in {engine}, importing...", "yellow"))
    load = subprocess.run([engine, "load", "-i", IMAGE_ARCHIVE], capture_output=True, text=True)
    if load.returncode != 0:
        print(_color(f"Failed to import GeneMarkS container: {load.stderr.strip()}"))
        return False
    print(_color("GeneMarkS container imported successfully", "green"))
    return True


def _run_container(engine, input_path, output_dir):
    fasta_dir = os.path.dirname(input_path)
    fasta = os.path.basename(input_path)
    image = IMAGE_NAME if engine == "podman" else f"localhost/{IMAGE_NAME}:latest"
    vol = ":Z" if _selinux_enforcing() else ""

    with tempfile.TemporaryDirectory(prefix="pelican_gmkey_") as tmp:
        key = os.path.join(tmp, "gm_key.gz")
        with gzip.open(key, "wb") as fh:
            fh.write(_read_key(find_gm_key()))
        os.chmod(tmp, 0o755)
        os.chmod(key, 0o644)

        cmd = [engine, "run", "--rm", "--entrypoint", ""]
        if engine == "docker":
            # Avoid root-owned output files and give the container a writable HOME
            cmd += ["--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp"]
        cmd += [
            "-v", f"{os.path.abspath(fasta_dir)}:/data{vol}",
            "-v", f"{os.path.abspath(output_dir)}:/output{vol}",
            "-v", f"{key}:/gm_key.gz{vol}",
            "-w", "/opt/gmsuite" if engine == "podman" else "/tmp",
            image,
            "bash", "-c",
            (
                f"exec > /output/{LOG_NAME} 2>&1 && "
                f"export LC_ALL=C && "  # Força idioma inglês
                f"echo 'Starting GeneMark processing...' && "
                f"echo 'Input file: {fasta}' && "
                f"cp /data/{fasta} . && "
                f"echo 'FASTA file copied successfully' && "
                f"gunzip -c /gm_key.gz > ~/.gm_key && "
                f"echo 'License key extracted successfully' && "
                f"echo 'Running GeneMark...' && "
                f"/opt/gmsuite/gmsn.pl --phage {fasta} && "
                f"echo 'GeneMark execution completed' && "
                f"cp {fasta}.lst /output/ && "
                f"echo 'Output file copied to /output/' && "
                f"echo 'Process finished successfully'"
            )
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"[{engine.upper()} returncode]: {result.returncode}")
    if result.stdout:
        print(f"[{engine.upper()} stdout]: {result.stdout}")
    if result.stderr:
        print(f"[{engine.upper()} stderr]: {result.stderr}")
    return result.returncode


###############
## Public API ##
###############

def select_backend(verbose=True):
    """Choose the GeneMarkS backend (cached). Returns None if none is usable."""
    global _selected_backend
    if _selected_backend:
        return _selected_backend

    requested = os.environ.get("PELICAN_GENEMARK_BACKEND", "auto").strip().lower()
    candidates = BACKENDS if requested in ("", "auto") else (requested,)
    if requested not in ("", "auto") and requested not in BACKENDS:
        print(_color(f"Unknown PELICAN_GENEMARK_BACKEND='{requested}'. Use: auto, native, podman or docker"))
        return None

    for backend in candidates:
        ok, reason = native_available() if backend == "native" else container_available(backend)
        if ok:
            _selected_backend = backend
            if verbose:
                print(_color(f"GeneMarkS backend: {backend}", "green"))
            return backend
        if verbose:
            print(_color(f"GeneMarkS backend '{backend}' unavailable: {reason}", "yellow"))
    return None


def prepare_genemark():
    """Select the backend and, for containers, make sure the image is loaded."""
    backend = select_backend()
    if backend is None:
        print(_color("No GeneMarkS backend available (native, podman or docker). "
                     "Run 'pelican-setup' for a diagnosis."))
        return None
    if backend in ("podman", "docker") and not ensure_image(backend):
        return None
    return backend


def run_podman(input_path, output_dir):
    """Run GeneMarkS (gmsn.pl --phage) on input_path; writes <fasta>.lst to output_dir."""
    if not os.path.isfile(input_path):
        print(f"FASTA file not found: {input_path}")
        return

    os.makedirs(output_dir, exist_ok=True)  # Garante que o diretório de saída existe
    input_path = os.path.abspath(input_path)
    output_lst = os.path.join(output_dir, os.path.basename(input_path) + ".lst")

    # Remove old files
    if os.path.exists(output_lst):
        os.remove(output_lst)

    backend = select_backend(verbose=False) or prepare_genemark()
    if backend is None:
        return
    if backend == "native":
        returncode = _run_native(input_path, output_dir)
    else:
        returncode = _run_container(backend, input_path, output_dir)

    log_path = os.path.join(output_dir, LOG_NAME)
    if license_expired(log_path):
        explain_license_error()
    elif returncode != 0:
        print(_color(f"GeneMarkS failed (backend: {backend}). See log: {log_path}"))
    return returncode
