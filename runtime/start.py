#!/opt/mamba/envs/latch-system/bin/python

import json
import os
import re
import sys
import urllib.request
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

sdk_token = (
    (latch_p / "token").read_text().strip() if (latch_p / "token").exists() else ""
)
skills_dir = Path("/opt/latch/plots-faas/.claude/skills")
skills_dir.mkdir(parents=True, exist_ok=True)

skills_branch = "main"
if sdk_token is not None:
    try:
        gql_body = json.dumps({
            "query": """
                query GetNotebookMetadata($podId: BigInt!) {
                    podInfo(id: $podId) {
                        plotNotebook { metadata }
                    }
                }
            """,
            "variables": {"podId": pod_id},
        }).encode()

        req = urllib.request.Request(
            f"https://vacuole.{domain}/graphql",
            data=gql_body,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Latch-SDK-Token {sdk_token}",
            },
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            metadata_str = (
                json
                .loads(resp.read())
                .get("data", {})
                .get("podInfo", {})
                .get("plotNotebook", {})
                .get("metadata")
            )
        if metadata_str is not None:
            skills_branch = json.loads(metadata_str).get("skillsBranch", "main")
    except Exception as e:
        print(f"failed to fetch notebook metadata: {e}", file=sys.stderr)

# todo: surface an error to the user if skills repo fails to pull or clone
latch_skills_dest = skills_dir / "latch-skills"
if latch_skills_dest.exists():
    ret = os.system(
        f"git -C {latch_skills_dest} fetch --depth 1 origin {skills_branch} && git -C {latch_skills_dest} checkout FETCH_HEAD"
    )
    if ret == 0:
        print(f"updated latch-skills to {skills_branch}")
    else:
        print(f"failed to update latch-skills to {skills_branch}", file=sys.stderr)
else:
    ret = os.system(
        f"git clone --depth 1 --branch {skills_branch} https://github.com/latchbio/latch-skills.git {latch_skills_dest}"
    )
    if ret == 0:
        print(f"cloned public latch-skills ({skills_branch}) -> {latch_skills_dest}")
    else:
        print("failed to clone public latch-skills", file=sys.stderr)

if latch_skills_dest.exists():
    for link in skills_dir.iterdir():
        if link.is_symlink() and str(latch_skills_dest) in str(link.readlink()):
            link.unlink()

    for sub in sorted(latch_skills_dest.iterdir()):
        if sub.is_dir() and (sub / "SKILL.md").exists():
            link = skills_dir / sub.name
            if not link.exists():
                link.symlink_to(sub)
                print(f"linked latch skill: {sub.name} -> {link}")

if sdk_token:
    try:
        gql_body = json.dumps({
            "query": """
                query AgentSkillRepos {
                    accountInfoCurrent {
                        id
                    }
                }
            """
        }).encode()

        req = urllib.request.Request(
            f"https://vacuole.{domain}/graphql",
            data=gql_body,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Latch-SDK-Token {sdk_token}",
            },
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            account_id = (
                json
                .loads(resp.read())
                .get("data", {})
                .get("accountInfoCurrent", {})
                .get("id")
            )

        if account_id is None:
            raise RuntimeError("could not resolve account")

        gql_body = json.dumps({
            "query": """
                query AgentSkillReposForAccount($accountId: BigInt!) {
                    agentSkillReposForAccount(argAccountId: $accountId)
                }
            """,
            "variables": {"accountId": account_id},
        }).encode()

        req = urllib.request.Request(
            f"https://vacuole.{domain}/graphql",
            data=gql_body,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Latch-SDK-Token {sdk_token}",
            },
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            skill_data = (
                json.loads(resp.read()).get("data", {}).get("agentSkillReposForAccount")
            )

        if skill_data is not None:
            pat = skill_data.get("pat")
            repos = skill_data.get("repos", [])

            if pat and repos:
                username = pat["username"]
                token = pat["token"]
                skills_dir.mkdir(parents=True, exist_ok=True)

                seen_dirs: set[str] = set()
                for repo in repos:
                    repo_url: str = repo["repoUrl"]
                    authed_url = repo_url.replace(
                        "https://", f"https://{username}:{token}@"
                    )

                    parts = repo_url.rstrip("/").removesuffix(".git").split("/")
                    repo_name = parts[-1] if parts else repo["displayName"]

                    if repo_name in seen_dirs:
                        print(
                            f"skill repo conflict: {repo_name} already cloned, skipping {repo_url}",
                            file=sys.stderr,
                        )
                        continue
                    seen_dirs.add(repo_name)

                    dest = skills_dir / repo_name
                    if not dest.exists():
                        ret = os.system(f"git clone --depth 1 {authed_url} {dest}")
                        if ret == 0:
                            print(f"cloned skill repo: {repo_url} -> {dest}")
                        else:
                            print(
                                f"failed to clone skill repo: {repo_url}",
                                file=sys.stderr,
                            )
                            continue

                    if not (dest / "SKILL.md").exists():
                        for sub in sorted(dest.iterdir()):
                            if sub.is_dir() and (sub / "SKILL.md").exists():
                                link = skills_dir / sub.name
                                if sub.name in seen_dirs:
                                    print(
                                        f"skill repo conflict: {sub.name} already exists, skipping {sub}",
                                        file=sys.stderr,
                                    )
                                    continue
                                seen_dirs.add(sub.name)
                                if not link.exists():
                                    link.symlink_to(sub)
                                    print(
                                        f"linked monorepo skill: {sub.name} -> {link}"
                                    )
    except Exception as e:
        print(f"failed to fetch skill repos: {e}", file=sys.stderr)

os.execle(
    "/usr/bin/nice",
    "-n",
    "-2",
    "/opt/mamba/envs/plots-faas/bin/python",
    "-m",
    "runtime.main",
    env_vars,
)
