from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from hercules.core.model.HerculesModel import HerculesModel, posterior_moments
from hercules.core.model.gridsearch.HyperparameterGrid import HyperparameterGrid


def test_stage1_beta_is_gamma_times_mu() -> None:
    gamma = np.array([0.25, 0.8, 1.0])
    mu = np.array([0.4, -0.2, 0.05])
    slab_variance = np.array([0.3, 0.1, 0.02])

    beta, _ = posterior_moments(gamma, mu, slab_variance)

    np.testing.assert_allclose(beta, gamma * mu)


def test_stage1_var_beta_is_marginal_posterior_variance() -> None:
    gamma = np.array([0.25, 0.8, 1.0])
    mu = np.array([0.4, -0.2, 0.05])
    slab_variance = np.array([0.3, 0.1, 0.02])

    _, variance = posterior_moments(gamma, mu, slab_variance)
    expected = gamma * (slab_variance + mu**2) - (gamma * mu) ** 2

    np.testing.assert_allclose(variance, expected)


def test_var_prior_precision_enters_initial_variational_parameters() -> None:
    model = object.__new__(HerculesModel)
    model.shapes = {"22": (2, 1)}
    model.Nj = {"22": np.array([[1000.0], [1000.0]])}
    model.sigma_epsilon = np.array([1.0])
    model.tau_beta = {"22": np.array([[50.0], [25.0]])}
    model.pi = np.array([0.1])
    model.float_precision = "float64"
    model.order = "F"

    model.initialize_variational_parameters()

    np.testing.assert_allclose(model.var_tau["22"][:, 0], [1050.0, 1025.0])


def test_grid_candidate_keeps_pi_and_sigma_epsilon_fixed() -> None:
    model = object.__new__(HerculesModel)
    model.fix_params = {
        "pi": np.array([0.01, 0.1]),
        "sigma_epsilon": np.array([0.8, 0.9]),
    }
    model.pi = np.array([0.01, 0.1])
    model.sigma_epsilon = np.array([0.8, 0.9])
    model.var_gamma = {"22": np.array([[0.5, 0.6], [0.7, 0.8]])}

    model.update_pi()
    model.update_sigma_epsilon()

    np.testing.assert_array_equal(model.pi, [0.01, 0.1])
    np.testing.assert_array_equal(model.sigma_epsilon, [0.8, 0.9])


def test_default_published_grid_has_one_hundred_candidates() -> None:
    grid = HyperparameterGrid(
        n_snps=10000,
        pi_steps=10,
        sigma_epsilon_steps=10,
    )

    table = grid.to_table()

    assert len(table) == 100
    assert table["pi"].nunique() == 10
    assert table["sigma_epsilon"].nunique() == 10

