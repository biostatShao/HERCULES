"""Command-line interface for HERCULES."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence

from . import __version__
from .config import load_config
from .diagnostics import collect_diagnostics, mandatory_runtime_ready
from .exceptions import HerculesError
from .workflow import execute_workflow, plan_workflow


def build_parser(*, program: str = "hercules") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=program, description="HERCULES cross-ancestry PRS workflow")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="report dependency availability")
    doctor.add_argument("--json", action="store_true", help="emit machine-readable diagnostics")
    doctor.add_argument("--config", help="inspect tool paths from a HERCULES configuration")
    doctor.set_defaults(handler=_doctor)

    config_parser = subparsers.add_parser("config", help="configuration operations")
    config_subparsers = config_parser.add_subparsers(dest="config_command", required=True)
    validate = config_subparsers.add_parser("validate", help="validate a YAML configuration")
    validate.add_argument("config")
    validate.add_argument(
        "--no-check-paths", action="store_true", help="validate structure without checking input paths"
    )
    validate.set_defaults(handler=_validate)

    run = subparsers.add_parser("run", help="run the dependency-safe full workflow")
    _add_execution_arguments(run)
    run.set_defaults(handler=_run, target="ensemble")

    stage = subparsers.add_parser("stage", help="run one dependency-safe model stage")
    stage_subparsers = stage.add_subparsers(dest="stage_name", required=True)
    for stage_name in ("m1", "m2", "m3"):
        stage_parser = stage_subparsers.add_parser(stage_name)
        _add_execution_arguments(stage_parser)
        stage_parser.set_defaults(handler=_run, target=stage_name)

    ensemble = subparsers.add_parser("ensemble", help="run the validated ensemble procedure")
    _add_execution_arguments(ensemble)
    ensemble.set_defaults(handler=_run, target="ensemble")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except HerculesError as exc:
        parser.error(str(exc))
    return 2


def _add_execution_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true", help="validate and print the workflow plan")
    parser.add_argument("--threads", type=int, help="override execution.threads")
    parser.add_argument("--parallel-jobs", type=int, help="override execution.parallel_jobs")
    parser.add_argument("--seed", type=int, help="override execution.seed")


def _overrides(args: argparse.Namespace) -> dict[str, object]:
    candidates = {
        "execution.threads": getattr(args, "threads", None),
        "execution.parallel_jobs": getattr(args, "parallel_jobs", None),
        "execution.seed": getattr(args, "seed", None),
    }
    return {key: value for key, value in candidates.items() if value is not None}


def _validate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    warnings_found = config.validate(check_paths=not args.no_check_paths)
    print(f"Configuration is valid for trait {config.trait_name!r}.")
    for message in warnings_found:
        print(f"WARNING: {message}")
    return 0


def _run(args: argparse.Namespace) -> int:
    config = load_config(args.config, cli_overrides=_overrides(args))
    warnings_found = config.validate(check_paths=not args.dry_run)
    for message in warnings_found:
        print(f"WARNING: {message}", file=sys.stderr)
    plan = plan_workflow(config, args.target)
    if args.dry_run:
        print(json.dumps(plan.as_dict(), indent=2))
        return 0
    execute_workflow(config, args.target)
    return 0


def _doctor(args: argparse.Namespace) -> int:
    config = load_config(args.config) if args.config else None
    diagnostics = collect_diagnostics(config.tools if config is not None else None)
    if args.json:
        print(json.dumps([item.as_dict() for item in diagnostics], indent=2))
    else:
        for item in diagnostics:
            marker = "OK" if item.available else "MISSING"
            print(f"[{marker:7}] {item.category:20} {item.name}: {item.detail}")
        print(
            "\nDependency availability permits execution but does not by itself demonstrate "
            "old-versus-new scientific equivalence."
        )
    return 0 if mandatory_runtime_ready(diagnostics) else 1
