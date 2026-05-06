import numpy as np

import genulens


def test_python_binding():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)
    cfg.model.imf.alpha2 = -1.13449983242887
    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    assert isinstance(arr, np.ndarray)
    assert arr.shape[0] == 10
    assert "tE" in result.columns


def test_python_callable_likelihood():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)

    def my_like(event):
        return 1.0 if event.tE > 10 else 0.0

    result = genulens.simulate(cfg, likelihood=my_like)
    arr = result.to_numpy()
    assert arr.shape[0] == 10
    assert np.all(arr[:, result.columns.index("tE")] > 10.0)


def test_python_typed_observation_config():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.observation.tE_obs = 20.0
    cfg.observation.tE_err = 0.1

    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    tE = arr[:, result.columns.index("tE")]
    assert tE.size > 0
    assert np.all((tE > 19.0) & (tE < 21.0))


def test_python_typed_model_config():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.model.density.add_x = 0
    cfg.model.density.stellar_halo = 0
    cfg.model.kinematics.omega_p = 45.0
    cfg.model.nsd.enabled = 0
    cfg.model.bh_kick.kick_bh = 50.0

    assert genulens.simulate(cfg).to_numpy().shape[0] == 5


def test_python_ruc_alias():
    result = genulens.ruc(l=0.0, b=0.0, n_simu=3, seed=1)
    assert result.to_numpy().shape[0] == 3


def test_python_isochrone_grid_lookup():
    grid = genulens.IsochroneGrid.load_default_roman()
    assert grid.row_count == 14723
    assert grid.sequence_count == 41
    assert "F146mag" in grid.bands

    query = genulens.IsochroneQuery()
    query.component = "thin1"
    query.log_age = 8.0
    query.metallicity_mh = -0.5
    query.initial_mass_msun = 0.1000000015

    star = grid.lookup(query)
    assert star.component == "thin1"
    assert np.isclose(star.teff_k, 2886.02441)
    assert np.isclose(star.absolute_magnitudes["F146mag"], 9.251)


def test_python_stellar_population_lookup():
    population = genulens.StellarPopulationModel.load_default_roman()
    assert genulens.StellarPopulationModel.component_name(7) == "thick"
    assert genulens.StellarPopulationModel.component_index("NSD") == 9

    query = genulens.StellarPopulationQuery()
    query.component_index = 7
    query.initial_mass_msun = 0.1000000015

    star = population.lookup(query)
    assert star.component == "thick"
    assert np.isclose(star.log_age, 10.07918)
    assert np.isclose(star.metallicity_mh, -0.8)


def test_python_forward_source_generator():
    generator = genulens.ForwardSourceGenerator.load_default_roman()
    query = genulens.ForwardSourceQuery()
    query.component = "thin1"
    query.distance_pc = 8000.0
    query.min_initial_mass_msun = 0.1
    query.max_initial_mass_msun = 0.11

    source = generator.sample(query, seed=123)
    assert source.stellar.component == "thin1"
    assert 0.1 <= source.stellar.initial_mass_msun <= 0.11
    assert source.angular_radius_microarcsec > 0.0
    assert np.isclose(
        source.apparent_magnitudes["F146mag"],
        source.stellar.absolute_magnitudes["F146mag"] + genulens.distance_modulus(8000.0),
    )
