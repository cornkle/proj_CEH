from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "configs"
DEFAULT_CONFIG = CONFIG_DIR / "paths.example.yml"
LOCAL_CONFIG = CONFIG_DIR / "paths.local.yml"


@dataclass
class PathConfig:
    machine: str
    roots: dict[str, str]
    datasets: dict[str, str]
    raw: dict[str, Any]

    def __getitem__(self, key: str) -> Any:
        if key in self.datasets:
            return self.datasets[key]
        if key in self.roots:
            return self.roots[key]
        if key in self.raw:
            return self.raw[key]
        raise KeyError(key)

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default


def _resolve_mapping(mapping: dict[str, Any], context: dict[str, str]) -> dict[str, Any]:
    resolved: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, str):
            value = os.path.expandvars(value)
            try:
                value = value.format_map(context)
            except KeyError:
                pass
            resolved[key] = os.path.normpath(value)
        elif isinstance(value, dict):
            resolved[key] = _resolve_mapping(value, context)
        else:
            resolved[key] = value
    return resolved


def load_paths(config_file: str | Path | None = None, allow_example: bool = False) -> PathConfig:
    if config_file is not None:
        config_path = Path(config_file)
    elif LOCAL_CONFIG.exists():
        config_path = LOCAL_CONFIG
    elif allow_example:
        config_path = DEFAULT_CONFIG
    else:
        raise FileNotFoundError(
            "No local paths config found. Create configs/paths.local.yml from configs/paths.example.yml."
        )

    if not config_path.exists():
        raise FileNotFoundError(f"Paths config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    machine = config.get("machine", "unknown")
    raw_roots = config.get("roots", {}) or {}
    roots = _resolve_mapping(raw_roots, {})
    raw_datasets = config.get("datasets", {}) or {}
    datasets = _resolve_mapping(raw_datasets, {**roots})

    return PathConfig(machine=machine, roots=roots, datasets=datasets, raw=config)
