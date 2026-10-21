"""Versioned *read-only* paper-method manifest.

No imports from frozen paper packages and no shell execution. CI checks only
manifest structure and checked-in source entry points; physics capabilities
are never inferred merely from a successful metadata validation.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_MANIFEST = HERE / "protocol_registry.json"
SCIENTIFIC_STATES = {"PASS", "PARTIAL", "BLOCKED", "NOT_STARTED"}
ALLOWED_CLAIMS = {"MODEL_VALIDATED_ONLY", "METHOD_VALIDATED_ONLY", "DIAGNOSTIC", "SYNTHETIC_ONLY"}
REQUIRED = {"id", "paper", "directory", "science_status", "install", "verify",
            "entrypoint", "default_args", "protocols", "capabilities", "limits", "evidence"}

class InvalidRegistry(ValueError):
    """Registry cannot support safe continuation."""

def _relative_path(value: str) -> Path:
    if not isinstance(value, str) or not value or "\\" in value:
        raise InvalidRegistry("invalid repository-relative path")
    p = Path(value)
    if p.is_absolute() or ".." in p.parts or "." in p.parts:
        raise InvalidRegistry("absolute or parent traversal is forbidden")
    return p

def load_registry(path: Path | None = None) -> dict:
    file = Path(path) if path is not None else DEFAULT_MANIFEST
    with file.open("r", encoding="utf-8") as stream:
        return json.load(stream)

def validate_registry(registry: dict, repo_root: Path | None = None,
                      check_files: bool = True) -> dict:
    root = (Path(repo_root) if repo_root is not None else REPO_ROOT).resolve()
    if registry.get("schema_version") != "1.0.0":
        raise InvalidRegistry("schema v1.0.0 required")
    if registry.get("registry_kind") != "metadata_only_no_physics_dispatch":
        raise InvalidRegistry("cannot claim a numerical engine")
    if registry.get("scientific_policy") != "CI_PASS_NOT_FULL_PAPER_PASS":
        raise InvalidRegistry("science acceptance guard required")
    modules = registry.get("modules")
    if not isinstance(modules, list) or len(modules) < 1:
        raise InvalidRegistry("modules required")
    ids = set()
    protocols = set()
    for m in modules:
        if not REQUIRED <= m.keys():
            raise InvalidRegistry("module missing required keys")
        mid = m["id"]
        if not isinstance(mid, str) or not mid or mid in ids:
            raise InvalidRegistry("duplicate/invalid module id")
        ids.add(mid)
        if m["science_status"] not in SCIENTIFIC_STATES:
            raise InvalidRegistry("invalid scientific status")
        dirname = _relative_path(m["directory"])
        if dirname.parts[0] != "papers":
            raise InvalidRegistry("paper registry may only reference papers/")
        full = (root / dirname).resolve()
        if not full.is_relative_to(root):
            raise InvalidRegistry("directory escapes root")
        if check_files and not full.is_dir():
            raise InvalidRegistry("paper directory absent: " + str(dirname))
        for k in ("install", "verify", "entrypoint", "default_args"):
            command = m[k]
            if not isinstance(command, list) or not all(isinstance(x, str) and x for x in command):
                raise InvalidRegistry("invalid command token list: "+k)
            if k != "default_args" and not command:
                raise InvalidRegistry("command may not be empty: "+k)
            if any(tok in command for tok in ("bash", "sh", "-c", "curl", "wget")):
                raise InvalidRegistry("unsafe manifest command")
        for command in (m["verify"], m["entrypoint"]):
            for tok in command:
                if tok.endswith(".py"):
                    p = _relative_path(tok)
                    if check_files and not (full/p).is_file():
                        raise InvalidRegistry(f"missing executable {dirname}/{tok}")
        doc = _relative_path(m["evidence"])
        if check_files and not (full/doc).is_file():
            raise InvalidRegistry("missing acceptance evidence: "+str(full/doc))
        if not isinstance(m["protocols"], list) or not m["protocols"]:
            raise InvalidRegistry("explicit protocols required")
        for proto in m["protocols"]:
            if not {"id", "kind", "status", "description"} <= proto.keys():
                raise InvalidRegistry("incomplete protocol")
            pid = proto["id"]
            if pid in protocols or not isinstance(pid, str) or "." not in pid:
                raise InvalidRegistry("protocols must be globally unique stable ids")
            protocols.add(pid)
            if proto["status"] not in ALLOWED_CLAIMS:
                raise InvalidRegistry("unsupported scientific claim")
            if not proto["description"]:
                raise InvalidRegistry("protocol scientific caveat absent")
        for key in ("capabilities", "limits"):
            if not isinstance(m[key], list) or not m[key]:
                raise InvalidRegistry("capability and limitation documentation required")
    for external in registry.get("external_dependencies", []):
        if "integration_status" not in external or not external.get("notes"):
            raise InvalidRegistry("external integration limits mandatory")
    return {"module_count":len(ids), "protocol_count":len(protocols),
            "source_file_check":bool(check_files),
            "status":"METADATA_VALID_ONLY",
            "note":"No paper test, DFT calculation or response integration was run"}

def get_module(module_id: str, registry: dict | None = None) -> dict:
    src = registry if registry is not None else load_registry()
    for module in src["modules"]:
        if module["id"] == module_id:
            return module
    raise KeyError(module_id)

def plan_module(module_id: str, protocol_id: str,
                registry: dict | None = None) -> dict:
    """Return a proposed *human reviewed* run; does not execute anything."""
    src = registry if registry is not None else load_registry()
    module = get_module(module_id, src)
    protocol = next((p for p in module["protocols"] if p["id"] == protocol_id), None)
    if protocol is None:
        raise InvalidRegistry("protocol not available for requested module")
    return {"module":module_id, "protocol":protocol,
            "cwd":module["directory"], "install_argv":module["install"],
            "verify_argv":module["verify"],
            "entrypoint_argv":module["entrypoint"]+module["default_args"],
            "execution":"NOT_EXECUTED",
            "validation":"REQUIRES_PAPER_TESTS_AND_PROTOCOL_CHECK",
            "science_status":module["science_status"],
            "limits":module["limits"]}
