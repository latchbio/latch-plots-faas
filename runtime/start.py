#!/opt/mamba/envs/latch-system/bin/python

import os
import re
import sys
from io import TextIOWrapper
from pathlib import Path

assert isinstance(sys.stdout, TextIOWrapper)
assert isinstance(sys.stderr, TextIOWrapper)

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

latch_p = Path("/root/.latch")

nucleus_url = (latch_p / "nucleus-url").read_text()
domain = ".".join(nucleus_url.split(".")[-2:])

pod_id = int((latch_p / "id").read_text())
plots_faas_version = Path("/opt/latch/plots_faas_version").read_text(encoding="utf-8")
dd_version = f"{pod_id}_{plots_faas_version}"

host_ip_pattern = re.compile(rb"HOST_IP=([^\0]+)")
proc_1_env_bytes = Path("/proc/1/environ").read_bytes()
host_ip_match = host_ip_pattern.search(proc_1_env_bytes)
host_ip = host_ip_match.group(1).decode() if host_ip_match else ""

otel_exporter_otlp_endpoint = f"http://{host_ip}:4317"

env_vars = {
    "DD_SERVICE": "plots-faas",
    "DD_ENV": "prod" if domain == "latch.bio" else "dev",
    "DD_VERSION": dd_version,
    "OTEL_EXPORTER_OTLP_ENDPOINT": otel_exporter_otlp_endpoint,
    "auto_reload": "false",
    "logging_mode": "console",
    "domain": domain,
    "auth_audience": "",
    "auth_self_signed_jwk": "",
    "CONDA_DEFAULT_ENV": "plots-faas",
    "CONDA_PREFIX": "/opt/mamba/envs/plots-faas",
    "PATH": "/opt/mamba/envs/plots-faas/bin:" + os.environ["PATH"],
    "PYTHON_GIL": "1",
}

os.system(
    "git -C /opt/latch/plots-faas remote add forgejo-mirror https://git.latch.bio/LatchBio/latch-plots-faas.git 2>/dev/null"
)

os.system(
    "git -C /opt/latch/plots-faas fetch --no-tags --depth 1 origin rteqs/start_script_optimization "
    "&& git -C /opt/latch/plots-faas reset --hard origin/rteqs/start_script_optimization "
    "|| (git -C /opt/latch/plots-faas fetch --no-tags --depth 1 forgejo-mirror main "
    "&& git -C /opt/latch/plots-faas reset --hard forgejo-mirror/main)"
)

os.system(
    "git -C /opt/latch/plots-faas submodule update --init --depth 1 --single-branch"
)
os.system("git -C /opt/latch/plots-faas rev-parse HEAD > /opt/latch/plots_faas_version")

os.chdir("/opt/latch/plots-faas")

os.system(
    "/opt/mamba/envs/plots-faas/bin/pip install --upgrade --upgrade-strategy only-if-needed latch"
)

os.execle(
    "/usr/bin/nice",
    "-n",
    "-2",
    "/opt/mamba/envs/plots-faas/bin/python",
    "-m",
    "runtime.main",
    env_vars,
)
