from __future__ import annotations


def test_public_import_is_lightweight() -> None:
    import hercules

    assert hercules.__version__ == "0.1.0.dev0"
    assert hercules.get_stage("m1").output_prefix == "HERCULES_M1"


def test_public_scientific_components_resolve() -> None:
    from hercules.inference import HerculesGridSearch, HerculesModel, vi_bayes_paper

    assert HerculesModel.__module__ == "hercules.core.model.HerculesModel"
    assert (
        HerculesGridSearch.__module__
        == "hercules.core.model.gridsearch.HerculesGridSearch"
    )
    assert callable(vi_bayes_paper)

