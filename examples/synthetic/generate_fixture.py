"""Generate the public deterministic two-ancestry HERCULES smoke fixture."""

from __future__ import annotations

import argparse
import copy
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


SAMPLES = 160
VARIANTS = 64


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "data")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)

    for ancestry, offset in (("base", 0), ("target", 1)):
        vcf = root / f"{ancestry}.vcf"
        dosage_matrix = _write_vcf(vcf, ancestry=ancestry, offset=offset)
        prefix = root / "genotypes" / ancestry / ancestry
        prefix.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                args.plink2,
                "--vcf",
                str(vcf),
                "--make-bed",
                "--threads",
                "1",
                "--out",
                str(prefix),
            ],
            check=True,
        )
        _write_validation_files(root, ancestry, offset, dosage_matrix)
        if ancestry == "target":
            for cohort in ("validation", "test"):
                cohort_prefix = root / "genotypes" / f"target_{cohort}" / f"target_{cohort}"
                cohort_prefix.parent.mkdir(parents=True, exist_ok=True)
                subprocess.run(
                    [
                        args.plink,
                        "--bfile",
                        str(prefix),
                        "--keep",
                        str(root / f"target.stage3.{cohort}.keep"),
                        "--make-bed",
                        "--out",
                        str(cohort_prefix),
                    ],
                    check=True,
                )
        _write_sumstats(root, ancestry, offset)

        from magenpy.GWADataLoader import GWADataLoader

        loader = GWADataLoader(
            bed_files=str(prefix),
            backend="plink",
            temp_dir=str(root / "tmp" / ancestry),
            threads=1,
            verbose=False,
        )
        loader.compute_ld(
            "windowed",
            output_dir=str(root / "ld" / ancestry),
            dtype="float32",
            window_size=VARIANTS,
        )
        loader.cleanup()

    _write_annotations(root)
    _write_configs(
        root,
        _resolve_executable(args.plink),
        _resolve_executable(args.plink2),
        _resolve_executable(args.rscript),
    )


def _resolve_executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is not None:
        return str(Path(resolved).resolve())
    path = Path(value).expanduser()
    if path.is_file():
        return str(path.resolve())
    raise FileNotFoundError(f"Executable not found: {value}")


def _write_vcf(path: Path, *, ancestry: str, offset: int) -> np.ndarray:
    samples = [f"{ancestry}_{index + 1:02d}" for index in range(SAMPLES)]
    rng = np.random.default_rng(7209 + offset)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=22,length=50818468>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]),
    ]
    genotype_codes = ("0/0", "0/1", "1/1")
    dosage_matrix = np.empty((VARIANTS, SAMPLES), dtype=np.int8)
    for variant in range(VARIANTS):
        allele_frequency = 0.10 + 0.035 * (variant % 10) + 0.01 * offset
        dosage = rng.binomial(2, min(allele_frequency, 0.48), size=SAMPLES)
        if np.all(dosage == dosage[0]):
            dosage[0], dosage[1] = 0, 1
        dosage_matrix[variant] = dosage
        genotypes = [genotype_codes[int(value)] for value in dosage]
        lines.append(
            "\t".join(
                [
                    "22",
                    str(1_000_000 + variant * 10_000),
                    f"rsHERC{variant + 1:03d}",
                    "A" if variant % 2 == 0 else "C",
                    "G" if variant % 2 == 0 else "T",
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    *genotypes,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dosage_matrix


def _write_validation_files(
    root: Path, ancestry: str, offset: int, dosage_matrix: np.ndarray
) -> None:
    weights = np.array(
        [
            ((-1.0) ** variant)
            * (0.006 + 0.001 * (variant % 8) + 0.0005 * offset)
            for variant in range(VARIANTS)
        ]
    )
    raw_signal = dosage_matrix.T @ weights
    genetic_signal = (raw_signal - raw_signal.mean()) / raw_signal.std(ddof=0)
    binary_latent = genetic_signal + 0.25 * np.sin(np.arange(SAMPLES) / 3.0)
    binary_threshold = float(np.median(binary_latent))
    rows = []
    quantitative_rows = []
    binary_rows = []
    keep_rows = []
    for index in range(SAMPLES):
        iid = f"{ancestry}_{index + 1:02d}"
        age = 35 + index
        sex = index % 2
        quantitative = (
            1.2 * genetic_signal[index]
            + 0.02 * age
            + 0.2 * sex
            + 0.05 * offset
            + 0.15 * math.sin(index / 3.0)
        )
        binary = int(binary_latent[index] > binary_threshold)
        rows.append((iid, quantitative, binary, age, sex))
        # PLINK 2 assigns FID=0 when importing the VCF sample IDs below.
        quantitative_rows.append(("0", iid, quantitative))
        binary_rows.append(("0", iid, binary))
        keep_rows.append(("0", iid))
    pd.DataFrame(rows, columns=["IID", "phenotype", "binary", "age", "sex"]).to_csv(
        root / f"{ancestry}.phenotype.tsv", sep="\t", index=False
    )
    pd.DataFrame(quantitative_rows).to_csv(
        root / f"{ancestry}.validation.quantitative.pheno",
        sep="\t",
        index=False,
        header=False,
    )
    pd.DataFrame(binary_rows).to_csv(
        root / f"{ancestry}.validation.binary.pheno",
        sep="\t",
        index=False,
        header=False,
    )
    pd.DataFrame(keep_rows).to_csv(
        root / f"{ancestry}.validation.keep", sep="\t", index=False, header=False
    )
    if ancestry == "target":
        phenotype_table = pd.DataFrame(
            rows, columns=["IID", "phenotype", "binary", "age", "sex"]
        )
        for cohort, indices in (
            ("validation", range(0, SAMPLES // 2)),
            ("test", range(SAMPLES // 2, SAMPLES)),
        ):
            selected = phenotype_table.iloc[list(indices)].copy()
            selected.to_csv(
                root / f"target.stage3.{cohort}.tsv", sep="\t", index=False
            )
            pd.DataFrame(
                [("0", iid) for iid in selected["IID"]]
            ).to_csv(
                root / f"target.stage3.{cohort}.keep",
                sep="\t",
                index=False,
                header=False,
            )
            if cohort == "validation":
                for column, label in (("phenotype", "quantitative"), ("binary", "binary")):
                    pd.DataFrame(
                        [("0", row.IID, getattr(row, column)) for row in selected.itertuples()]
                    ).to_csv(
                        root / f"target.stage3.validation.{label}.pheno",
                        sep="\t",
                        index=False,
                        header=False,
                    )


def _write_sumstats(root: Path, ancestry: str, offset: int) -> None:
    rows = []
    for variant in range(VARIANTS):
        beta = ((-1.0) ** variant) * (
            0.006 + 0.001 * (variant % 8) + 0.0005 * offset
        )
        se = 0.02 + 0.001 * (variant % 4)
        p_value = math.erfc(abs(beta / se) / math.sqrt(2.0))
        prior_variance = 0.02 + 0.0005238095238095238 * variant + 0.001 * offset
        rows.append(
            (
                22,
                f"rsHERC{variant + 1:03d}",
                1_000_000 + variant * 10_000,
                "G" if variant % 2 == 0 else "T",
                "A" if variant % 2 == 0 else "C",
                50_000 + 2_500 * offset,
                0.15 + 0.02 * (variant % 5),
                beta,
                se,
                p_value,
                prior_variance,
            )
        )
    pd.DataFrame(
        rows,
        columns=[
            "CHR",
            "SNP",
            "POS",
            "A1",
            "A2",
            "N",
            "AF1",
            "BETA",
            "SE",
            "P",
            "var_prior",
        ],
    ).to_csv(root / f"{ancestry}.fastGWA.tsv", sep="\t", index=False)


def _write_annotations(root: Path) -> None:
    table = pd.DataFrame(
        {
            "SNP": [f"rsHERC{index + 1:03d}" for index in range(VARIANTS)],
            "baseline": np.linspace(0.02, 0.053, VARIANTS),
            "coding": [int(index % 3 == 0) for index in range(VARIANTS)],
            "conserved": [int(index % 4 < 2) for index in range(VARIANTS)],
        }
    )
    table.to_csv(root / "functional_annotation.tsv", sep="\t", index=False)
    table.loc[:, ["SNP", "baseline"]].to_csv(
        root / "per_snp_heritability.tsv", sep="\t", index=False
    )


def _write_configs(root: Path, plink: str, plink2: str, rscript: str) -> None:
    common = {
        "trait_name": "synthetic_quantitative",
        "chromosomes": ["22"],
        "base_ancestry": "BASE",
        "target_ancestry": "TARGET",
        "inputs": {
            "summary_statistics": {
                "base_path": str(root / "base.fastGWA.tsv"),
                "target_path": str(root / "target.fastGWA.tsv"),
                "base_columns": {},
                "target_columns": {},
            },
            "functional_annotation": str(root / "functional_annotation.tsv"),
            "per_snp_heritability": str(root / "per_snp_heritability.tsv"),
            "ld_reference": {
                "base": str(root / "ld" / "base" / "chr_22"),
                "target": str(root / "ld" / "target" / "chr_22"),
            },
            "genotype_prefixes": {
                "base_validation": str(root / "genotypes" / "base" / "base"),
            },
            "validation_genotype": str(
                root / "genotypes" / "target_validation" / "target_validation"
            ),
            "phenotype_file": "",
            "target_validation_genotype": str(
                root / "genotypes" / "target_validation" / "target_validation"
            ),
            "target_validation_phenotype": str(root / "target.stage3.validation.tsv"),
            "target_test_genotype": str(
                root / "genotypes" / "target_test" / "target_test"
            ),
            "target_test_phenotype": str(root / "target.stage3.test.tsv"),
            "phenotype_column": "phenotype",
            "covariates": ["age", "sex"],
            "trait_type": "quantitative",
        },
        "output_dir": str(root / "results" / "quantitative"),
        "temporary_dir": str(root / "results" / "quantitative" / "tmp"),
        "tools": {"plink": plink, "plink2": plink2, "rscript": rscript},
        "execution": {"threads": 1, "parallel_jobs": 1, "seed": 7209},
        "m1": {
            "hyperparameter_search": "grid",
            "pi_steps": 10,
            "sigma_epsilon_steps": 10,
            "sumstats_format": "fastgwa",
            "validation_keep": str(root / "target.stage3.validation.keep"),
            "validation_phenotype": str(
                root / "target.stage3.validation.quantitative.pheno"
            ),
        },
        "m2": {
            "hyperparameter_search": "grid",
            "pi_steps": 10,
            "sigma_epsilon_steps": 10,
            "sumstats_format": "fastgwa",
            "validation_genotype": str(root / "genotypes" / "base" / "base"),
            "validation_keep": str(root / "base.validation.keep"),
            "validation_phenotype": str(root / "base.validation.quantitative.pheno"),
        },
        "m3": {
            "model": "directional_pairwise_uniform_lambda",
            "lambda_prior": "uniform_0_1",
            "max_iter": 1000,
            "tol": 1e-6,
            "quadrature_points": 32,
        },
        "ensemble": {
            "quantitative_learners": ["lasso", "ridge", "neural_network"],
            "binary_learners": ["lasso", "neural_network"],
            "binary_method": "method.AUC",
        },
        "checkpoint": {"enabled": True, "resume": True},
        "logging": {"level": "INFO", "file": ""},
    }
    (root / "hercules.quantitative.yaml").write_text(
        yaml.safe_dump(common, sort_keys=False), encoding="utf-8"
    )
    binary = copy.deepcopy(common)
    binary["trait_name"] = "synthetic_binary"
    binary["inputs"] = dict(common["inputs"])
    binary["inputs"]["phenotype_column"] = "binary"
    binary["inputs"]["covariates"] = []
    binary["inputs"]["trait_type"] = "binary"
    binary["output_dir"] = str(root / "results" / "binary")
    binary["temporary_dir"] = str(root / "results" / "binary" / "tmp")
    binary["m1"]["validation_phenotype"] = str(
        root / "target.stage3.validation.binary.pheno"
    )
    binary["m2"]["validation_phenotype"] = str(
        root / "base.validation.binary.pheno"
    )
    (root / "hercules.binary.yaml").write_text(
        yaml.safe_dump(binary, sort_keys=False), encoding="utf-8"
    )


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    main()
