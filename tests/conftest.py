from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest


@pytest.fixture
def config_mapping(tmp_path: Path) -> dict[str, Any]:
    input_dir = tmp_path / "inputs with spaces"
    input_dir.mkdir()
    base = input_dir / "base.tsv"
    target = input_dir / "target.tsv"
    phenotype = input_dir / "phenotype.tsv"
    annotation = input_dir / "annotation.tsv"
    per_snp = input_dir / "per-snp.tsv"
    validation_pheno = input_dir / "validation.pheno"
    validation_keep = input_dir / "validation.keep"
    bed_prefix = input_dir / "validation cohort"
    target_validation_prefix = input_dir / "target validation cohort"
    target_test_prefix = input_dir / "target test cohort"
    target_validation_phenotype = input_dir / "target validation.tsv"
    target_test_phenotype = input_dir / "target test.tsv"
    ld_base = input_dir / "ld base"
    ld_target = input_dir / "ld target"
    for path in (
        base,
        target,
        phenotype,
        annotation,
        per_snp,
        validation_pheno,
        validation_keep,
        target_validation_phenotype,
        target_test_phenotype,
    ):
        path.write_text("fixture\n", encoding="utf-8")
    fastgwa_fixture = (
        "CHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP\n"
        "1\trs1\t1000\tA\tG\t1000\t0.25\t0.01\t0.02\t0.617075\n"
    )
    base.write_text(fastgwa_fixture, encoding="utf-8")
    target.write_text(fastgwa_fixture, encoding="utf-8")
    for prefix in (bed_prefix, target_validation_prefix, target_test_prefix):
        for suffix in (".bed", ".bim", ".fam"):
            Path(f"{prefix}{suffix}").write_bytes(b"")
    ld_base.mkdir()
    ld_target.mkdir()
    tool_dir = tmp_path / "tools with spaces"
    tool_dir.mkdir()
    plink = tool_dir / "plink.exe"
    plink2 = tool_dir / "plink2.exe"
    rscript = tool_dir / "Rscript.exe"
    for path in (plink, plink2, rscript):
        path.write_bytes(b"")
        path.chmod(0o755)
    return {
        "trait_name": "height",
        "chromosomes": ["1", "2"],
        "base_ancestry": "EUR",
        "target_ancestry": "AFR",
        "inputs": {
            "summary_statistics": {
                "base_path": str(base),
                "target_path": str(target),
                "base_columns": {"snp": "SNP", "beta": "BETA"},
                "target_columns": {"snp": "SNP", "beta": "BETA"},
            },
            "functional_annotation": str(annotation),
            "per_snp_heritability": str(per_snp),
            "ld_reference": {"base": str(ld_base), "target": str(ld_target)},
            "genotype_prefixes": {
                "base_validation": str(bed_prefix),
                "target": str(bed_prefix),
            },
            "validation_genotype": str(bed_prefix),
            "phenotype_file": str(phenotype),
            "target_validation_genotype": str(target_validation_prefix),
            "target_validation_phenotype": str(target_validation_phenotype),
            "target_test_genotype": str(target_test_prefix),
            "target_test_phenotype": str(target_test_phenotype),
            "phenotype_column": "height",
            "trait_type": "quantitative",
        },
        "output_dir": str(tmp_path / "output with spaces"),
        "temporary_dir": str(tmp_path / "temporary with spaces"),
        "tools": {"plink": str(plink), "plink2": str(plink2), "rscript": str(rscript)},
        "m1": {
            "validation_phenotype": str(validation_pheno),
            "validation_keep": str(validation_keep),
        },
        "m2": {
            "validation_genotype": str(bed_prefix),
            "validation_phenotype": str(validation_pheno),
            "validation_keep": str(validation_keep),
        },
    }


@pytest.fixture
def config_path(tmp_path: Path, config_mapping: dict[str, Any]) -> Path:
    path = tmp_path / "hercules config.yaml"
    # JSON is valid YAML and keeps the fixture usable if PyYAML is unavailable.
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    return path
