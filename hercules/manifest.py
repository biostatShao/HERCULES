"""Configuration-aware run manifests and checkpoint metadata."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .paths import safe_path_component
from .stages import StageSpec


def configuration_hash(config: Mapping[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class RunManifest:
    project: str
    version: str
    trait: str
    configuration_hash: str
    created_at: str
    stages: tuple[str, ...]

    @classmethod
    def create(
        cls, *, version: str, trait: str, config: Mapping[str, Any], stages: tuple[str, ...]
    ) -> RunManifest:
        return cls(
            project="HERCULES",
            version=version,
            trait=trait,
            configuration_hash=configuration_hash(config),
            created_at=datetime.now(timezone.utc).isoformat(),
            stages=stages,
        )

    def write(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(asdict(self), indent=2) + "\n", encoding="utf-8")
        return destination


def checkpoint_path(
    output_dir: str | Path,
    *,
    stage: StageSpec,
    trait: str,
    chromosome: str | int | None,
) -> Path:
    root = Path(output_dir).resolve()
    safe_trait = safe_path_component(trait, "trait")
    chrom = "all" if chromosome is None else f"chr{safe_path_component(str(chromosome), 'chromosome')}"
    candidate = root / "checkpoints" / safe_trait / chrom / stage.checkpoint_name
    if not candidate.resolve().is_relative_to(root):
        raise ValueError(f"Checkpoint path escapes output directory: {candidate}")
    return candidate


def checkpoint_matches(path: str | Path, expected_hash: str) -> bool:
    marker = Path(path)
    if not marker.is_file():
        return False
    try:
        data = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    if not isinstance(data, dict):
        return False
    return data.get("configuration_hash") == expected_hash and data.get("status") == "complete"
