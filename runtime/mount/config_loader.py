import importlib.util
import sys
from pathlib import Path
from textwrap import dedent

config_path = Path("/opt/latch/latch-plots-faas/runtime/mount/agent_config")

def load_config_module(module_name: str):
    config_file = config_path / f"{module_name}.py"

    if not config_file.exists():
        raise FileNotFoundError(f"Config module not found: {config_file}")

    spec = importlib.util.spec_from_file_location(f"agent_config.{module_name}", config_file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[f"agent_config.{module_name}"] = module
    spec.loader.exec_module(module)

    return module

def get_custom_tools(function_tool) -> list:
    try:
        config_file = config_path / "tools.py"

        if not config_file.exists():
            return []

        spec = importlib.util.spec_from_file_location("agent_config.tools", config_file)
        module = importlib.util.module_from_spec(spec)

        sys.modules["agents"] = type(sys)("agents")
        sys.modules["agents"].function_tool = function_tool

        sys.modules["agent_config.tools"] = module
        spec.loader.exec_module(module)

        return module.custom_tools
    except Exception as e:
        print(f"[config_loader] Warning: Could not load custom tools: {e}")
        return []

def build_full_instruction(initial_context: str) -> str:
    try:
        prompts = load_config_module("prompts")
    except Exception as e:
        print(f"[config_loader] Error: Could not load prompts module: {e}")
        raise

    loaded_docs = {}

    for doc_config in prompts.external_docs:
        name = doc_config["name"]
        path = Path(doc_config["path"])
        doc_type = doc_config["type"]

        try:
            if doc_type == "directory":
                content = ""
                if path.exists():
                    for file in path.rglob("*.mdx"):
                        doc_name = file.stem
                        file_content = file.read_text().strip()
                        content += f"<{name}_{doc_name}>\n{file_content}\n</{name}_{doc_name}>\n\n"
                loaded_docs[name] = content
            elif doc_type == "file":
                if path.exists():
                    loaded_docs[name] = path.read_text()
                else:
                    loaded_docs[name] = ""
        except Exception as e:
            print(f"[config_loader] Warning: Could not load {name} from {path}: {e}")
            loaded_docs[name] = ""

    instruction_parts = []

    for name, content in loaded_docs.items():
        if content:
            instruction_parts.append(f"<{name}>\n{content}\n</{name}>\n")

    instruction_parts.append(prompts.system_instruction)

    instruction_parts.append(f"\n<notebook_context>\n{initial_context}\n</notebook_context>\n")

    return dedent("\n".join(instruction_parts))
