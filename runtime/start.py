#!/opt/mamba/envs/latch-system/bin/python

import os
import re
from pathlib import Path

latch_p = Path("/root/.latch")

nucleus_url = (latch_p / "nucleus-url").read_text()
domain = ".".join(nucleus_url.split(".")[-2:])

pod_id = int((latch_p / "id").read_text())
plots_faas_version = Path("/opt/latch/plots_faas_version").read_text()
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
    "PATH": "/opt/latch/bin:/opt/mamba/envs/plots-faas/bin:/opt/mamba/condabin:/usr/local/bin:/usr/bin",
}

os.system("git -C /opt/latch/plots-faas pull origin main")
os.system("git -C /opt/latch/plots-faas rev-parse HEAD > /opt/latch/plots_faas_version")

os.chdir("/opt/latch/plots-faas")
os.execle(
    "/usr/bin/nice",
    "-n",
    "-2",
    "/opt/mamba/envs/plots-faas/bin/python",
    "-m",
    "runtime.main",
    env_vars,
)
