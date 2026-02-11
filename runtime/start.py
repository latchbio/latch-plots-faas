#!/opt/mamba/envs/latch-system/bin/python

import os
import re
from pathlib import Path

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
    "PYTHON_GIL": "0",
}

os.system(
    "git -C /opt/latch/plots-faas remote add forgejo-mirror https://git.latch.bio/LatchBio/latch-plots-faas.git"
)

os.system(
    "git -C /opt/latch/plots-faas pull origin main || git -C /opt/latch/plots-faas forgejo-mirror main"
)
os.system("git -C /opt/latch/plots-faas submodule update --init --remote")
os.system("git -C /opt/latch/plots-faas rev-parse HEAD > /opt/latch/plots_faas_version")

os.chdir("/opt/latch/plots-faas")

os.system("/opt/mamba/envs/plots-faas/bin/pip install --upgrade latch")

technology_docs = Path(
    "/opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs"
)
for sp_requirements in technology_docs.glob("*/requirements.txt"):
    os.system(
        f"/opt/mamba/envs/plots-faas/bin/pip install -c /opt/latch/requirements-plots.txt -r {sp_requirements}"
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
