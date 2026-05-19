import numpy as np
import pytest
from pathlib import Path

import genulens


def test_python_binding():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)
    cfg.model.imf.alpha2 = -1.13449983242887
    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    assert isinstance(arr, np.ndarray)
    assert arr.shape[0] == 10
    assert "t_E" in result.columns
    assert "pi_EN" in result.columns
    assert "mu_rel_N" in result.columns


def test_python_callable_likelihood():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)

    def my_like(event):
        return 1.0 if event.t_E > 10 else 0.0

    result = genulens.simulate(cfg, likelihood=my_like)
    arr = result.to_numpy()
    assert arr.shape[0] == 10
    assert np.all(arr[:, result.columns.index("t_E")] > 10.0)


def test_python_typed_observation_config():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.observation.tE_obs = 20.0
    cfg.observation.tE_err = 0.1

    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    tE = arr[:, result.columns.index("t_E")]
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
    assert grid.row_count == 41039
    assert grid.sequence_count == 113
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


def test_python_isochrone_grid_skips_discontinuous_segments():
    grid = genulens.IsochroneGrid.load_default_prime()
    query = genulens.IsochroneQuery()
    query.component = "bar"
    query.log_age = 9.903089987
    query.metallicity_mh = -0.5
    query.initial_mass_msun = 0.99758

    star = grid.lookup(query)
    assert star.teff_k < 4000.0
    assert star.absolute_magnitudes["Imag"] < -3.0


def test_python_stellar_population_lookup():
    population = genulens.StellarPopulationModel.load_default_roman()
    assert genulens.StellarPopulationModel.component_name(7) == "thick"
    assert genulens.StellarPopulationModel.component_index("NSD") == 9

    query = genulens.StellarPopulationQuery()
    query.component_index = 7
    query.initial_mass_msun = 0.1000000015

    star = population.lookup(query)
    assert star.component == "thick"
    assert np.isclose(star.log_age, 10.079181246047625)
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
    assert "F146mag" in source.stellar.absolute_magnitudes

    result = generator.sample_many(query, n_sources=4, seed=123)
    arr = result.to_numpy()
    assert arr.shape == (4, len(result.columns))
    assert "M_F146mag_S" in result.columns
    assert "F146mag" in result.bands


def test_python_simulation_forward_source_mode():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "roman"
    cfg.source.use_magnitude_selection = 0
    cfg.source.min_initial_mass_msun = 0.1
    cfg.source.max_initial_mass_msun = 0.2

    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    assert arr.shape[0] == 5
    assert "teff_S" in result.columns
    assert "M_F146mag_S" in result.columns
    teff = arr[:, result.columns.index("teff_S")]
    assert np.all(np.isfinite(teff))


def test_python_source_mode_convenience_methods():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.use_isochrone_source(
        i_min=12.0,
        i_max=21.0,
        band="Imag",
        photometry="prime",
        min_mass=0.1,
        max_mass=2.0,
    )
    cfg.use_genstars_extinction(dm_rc=14.5, ejk_rc=1.0)
    cfg.observation.IL_err = 0.0

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5
    assert "M_Imag_S" in result.columns


def test_python_simulation_forward_source_matches_source_selection():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "prime"
    cfg.source.band = "Ksmag_2mass"
    cfg.source.min_magnitude = 12.0
    cfg.source.max_magnitude = 25.0
    cfg.source.min_initial_mass_msun = 0.1
    cfg.source.max_initial_mass_msun = 1.0
    cfg.source.ak_rc = 0.2

    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    assert arr.shape[0] == 5
    assert "M_Ksmag_2mass_S" in result.columns
    assert np.all(np.isfinite(arr[:, result.columns.index("M_Ksmag_2mass_S")]))


def test_python_simulation_forward_source_h_band_manual_extinction():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "prime"
    cfg.source.band = "Hmag_2mass"
    cfg.source.min_magnitude = 12.0
    cfg.source.max_magnitude = 25.0
    cfg.source.min_initial_mass_msun = 0.1
    cfg.source.max_initial_mass_msun = 1.0
    cfg.source.dm_rc = 14.5
    cfg.source.ah_rc = 0.2

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5
    assert "M_Hmag_2mass_S" in result.columns


def test_python_simulation_forward_source_imag_auto_extinction_changes_source_prior():
    def base_config():
        cfg = genulens.Config(l=1.0, b=-3.9, n_simu=20, seed=1234)
        cfg.source.i_min = 12.0
        cfg.source.i_max = 21.0
        cfg.source.extinction_mode = "genstars"
        cfg.source.extinction_law = 1
        cfg.source.ejk_rc = 1.0
        cfg.source.dm_rc = 14.5
        cfg.observation.IL_err = 0.0
        return cfg

    legacy = genulens.simulate(base_config())

    cfg = base_config()
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "prime"
    cfg.source.band = "Imag"
    cfg.source.min_magnitude = 12.0
    cfg.source.max_magnitude = 21.0
    cfg.source.min_initial_mass_msun = 0.1
    cfg.source.max_initial_mass_msun = 2.0

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 20
    assert "M_Imag_S" in result.columns

    corner_columns = ["M_L", "D_L", "D_S", "mu_rel_N", "mu_rel_E"]
    legacy_arr = legacy.to_numpy()
    result_arr = result.to_numpy()
    legacy_idx = [legacy.columns.index(c) for c in corner_columns]
    result_idx = [result.columns.index(c) for c in corner_columns]
    assert not np.allclose(legacy_arr[:, legacy_idx], result_arr[:, result_idx])

    d_s = result_arr[:, result.columns.index("D_S")]
    abs_i = result_arr[:, result.columns.index("M_Imag_S")]
    ai_rc = 3.97
    dmean = 10 ** (0.2 * cfg.source.dm_rc) * 10
    hscale = cfg.source.dust_scale_height_pc / (
        abs(np.sin(np.deg2rad(cfg.b))) + 0.0001
    )
    ai0 = ai_rc / (1.0 - np.exp(-dmean / hscale))
    ai = ai0 * (1.0 - np.exp(-d_s / hscale))
    apparent_i = abs_i + 5.0 * np.log10(0.1 * d_s) + ai
    assert np.all((apparent_i >= 12.0) & (apparent_i <= 21.0))


def test_python_simulation_forward_source_h_band_requires_extinction():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "prime"
    cfg.source.band = "Hmag_2mass"
    cfg.source.min_magnitude = 12.0
    cfg.source.max_magnitude = 25.0

    with pytest.raises(RuntimeError, match="H band requires AHrc"):
        genulens.simulate(cfg)


def test_python_simulation_forward_source_absolute_selection():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.source.mode = "isochrone"
    cfg.source.photometry = "roman"
    cfg.source.band = "F146mag"
    cfg.source.min_magnitude = 3.0
    cfg.source.max_magnitude = 8.5
    cfg.source.apparent_magnitude = 0

    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    f146 = arr[:, result.columns.index("M_F146mag_S")]
    assert np.all((f146 >= 3.0) & (f146 <= 8.5))


def test_forward_source_isochrone_mixture_same_table_matches_default():
    table = "source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat"
    default = genulens.ForwardSourceGenerator.load_default_roman()
    mixed = genulens.ForwardSourceGenerator.load_mixture(table, table, 0.5)

    cut = genulens.MagnitudeSelection()
    cut.band = "F146mag"
    cut.min_magnitude = 3.0
    cut.max_magnitude = 8.5

    query = genulens.ForwardSourceQuery()
    query.component_index = 0
    query.distance_pc = 8000.0
    query.min_initial_mass_msun = 0.1
    query.max_initial_mass_msun = 1.0
    query.use_default_log_age = False
    query.log_age = 8.0
    query.use_default_metallicity = False
    query.metallicity_mh = 0.0
    query.magnitude_selections = [cut]

    assert np.isclose(
        mixed.selection_probability(query),
        default.selection_probability(query),
    )
    source = mixed.sample(query, seed=1234)
    assert np.isfinite(source.stellar.teff_k)


def test_python_default_isochrone_table_paths():
    spec = genulens.IsochroneLibrarySpec()
    spec.family = "parsec"
    spec.photometry = "roman"
    spec.abundance = "solar_scaled"
    assert (
        genulens.default_isochrone_table_path(spec)
        == "source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat"
    )

    spec.family = "mist"
    spec.abundance = "alpha_enhanced"
    spec.alpha_fe = 0.4
    assert (
        genulens.default_isochrone_table_path(spec)
        == "source_photometry/mist/v2.5/alpha_enhanced/afe_p0p40/normalized/all_roman_mist_alpha_afe_p0p40.dat"
    )

    spec.family = "parsec"
    spec.abundance = "alpha_enhanced"
    with pytest.raises(RuntimeError, match="PARSEC/CMD alpha-enhanced"):
        genulens.default_isochrone_table_path(spec)

    spec.family = "not_a_family"
    spec.abundance = "solar_scaled"
    with pytest.raises(RuntimeError, match="expected parsec or mist"):
        genulens.default_isochrone_table_path(spec)


def test_python_mist_alpha_table_covers_all_source_components():
    spec = genulens.IsochroneLibrarySpec()
    spec.family = "mist"
    spec.photometry = "roman"
    spec.abundance = "alpha_enhanced"
    spec.alpha_fe = 0.4
    grid = genulens.IsochroneGrid.load(Path(genulens.default_isochrone_table_path(spec)))

    component_names = [
        "thin1",
        "thin2",
        "thin3",
        "thin4",
        "thin5",
        "thin6",
        "thin7",
        "thick",
        "bar",
        "NSD",
        "halo",
    ]
    assert grid.sequence_count >= len(component_names)
    for component_index, component in enumerate(component_names):
        point = genulens.SourcePopulationPrior.points_for_component(component_index)[0]
        query = genulens.IsochroneQuery()
        query.component = component
        query.log_age = point.log_age
        query.metallicity_mh = point.metallicity_mh
        query.initial_mass_msun = 0.5
        star = grid.lookup(query)
        assert star.component_index == component_index
        assert np.isfinite(star.teff_k)
        assert "F146mag" in star.absolute_magnitudes


def test_python_simulation_forward_source_alpha_fraction_default_path_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.isochrone_model = "alpha_mixture"
    cfg.forward_source.alpha_enhanced_fraction = 0.5

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5


def test_python_simulation_forward_source_alpha_fraction_same_table_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.isochrone_model = "alpha_mixture"
    cfg.forward_source.alpha_enhanced_fraction = 0.5
    cfg.forward_source.alpha_enhanced_table_path = (
        "source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat"
    )

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5
    assert "teff_S" in result.columns


def test_python_simulation_forward_source_family_selection_same_table_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.isochrone_family = "parsec"
    cfg.forward_source.isochrone_abundance = "solar_scaled"

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5


def test_python_simulation_forward_source_mist_default_path_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.isochrone_family = "mist"
    cfg.forward_source.isochrone_abundance = "solar_scaled"

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5


def test_python_simulation_forward_source_component_mist_alpha_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.isochrone_model = "alpha_mixture"
    cfg.forward_source.secondary_isochrone_family = "mist"
    cfg.forward_source.secondary_isochrone_abundance = "alpha_enhanced"
    cfg.forward_source.secondary_isochrone_alpha_fe = 0.4
    cfg.forward_source.alpha_enhanced_fraction = 0.0
    cfg.forward_source.alpha_enhanced_components = [8, 9]
    cfg.forward_source.alpha_enhanced_component_fractions = [1.0, 1.0]

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5


def test_python_simulation_forward_source_component_alpha_fraction_same_table_runs():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.isochrone_model = "alpha_mixture"
    cfg.forward_source.alpha_enhanced_table_path = (
        "source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat"
    )
    cfg.forward_source.alpha_enhanced_fraction = 0.0
    cfg.forward_source.alpha_enhanced_components = [8, 9]
    cfg.forward_source.alpha_enhanced_component_fractions = [0.5, 1.0]

    result = genulens.simulate(cfg)
    assert result.to_numpy().shape[0] == 5


def test_python_simulation_forward_source_component_alpha_fraction_validates_lengths():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.isochrone_model = "alpha_mixture"
    cfg.forward_source.alpha_enhanced_table_path = (
        "source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat"
    )
    cfg.forward_source.alpha_enhanced_components = [8, 9]
    cfg.forward_source.alpha_enhanced_component_fractions = [0.5]

    with pytest.raises(RuntimeError, match="same length"):
        genulens.simulate(cfg)


def test_python_simulation_forward_source_apparent_requires_extinction():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "roman"
    cfg.forward_source.selection_bands = ["F146mag"]
    cfg.forward_source.selection_min_magnitudes = [17.0]
    cfg.forward_source.selection_max_magnitudes = [24.0]

    with pytest.raises(RuntimeError, match="requires auto extinction"):
        genulens.simulate(cfg)


def test_python_simulation_forward_source_zero_density_errors():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.forward_source.enabled = 1
    cfg.forward_source.photometry = "prime"
    cfg.source.extinction_mode = "genstars"
    cfg.source.extinction_law = 1
    cfg.source.ejk_rc = 1.0
    cfg.source.dm_rc = 14.5
    cfg.observation.IL_err = 0.0
    cfg.forward_source.selection_bands = ["Imag"]
    cfg.forward_source.selection_min_magnitudes = [-10.0]
    cfg.forward_source.selection_max_magnitudes = [-5.0]

    with pytest.raises(RuntimeError, match="zero source density"):
        genulens.simulate(cfg)


def test_python_simulation_cli_verbosity_columns():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    result = genulens.simulate(cfg)
    assert result.columns[:19] == [
        "wtj",
        "M_L",
        "D_L",
        "D_S",
        "t_E",
        "theta_E",
        "pi_E",
        "pi_EN",
        "pi_EE",
        "mu_rel",
        "mu_rel_N",
        "mu_rel_E",
        "mu_Sl",
        "mu_Sb",
        "I_L",
        "K_L",
        "iS",
        "iL",
        "fREM",
    ]
    arr = result.to_numpy()
    assert arr.shape == (5, 19)
    assert np.all(np.isfinite(arr[:, result.columns.index("mu_rel_N")]))


def test_python_simulation_verbosity_1_columns():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=5, seed=1234)
    cfg.sampling.verbosity = 1
    result = genulens.simulate(cfg)
    assert result.columns == [
        "wtj",
        "tE",
        "thetaE",
        "piE",
        "M_L",
        "D_S",
        "D_L",
        "mu_rel",
        "iS",
        "iL",
        "tau_s",
        "tau_l",
        "fREM",
    ]


def test_selection_probability_properties():
    generator = genulens.ForwardSourceGenerator.load_default_roman()

    query = genulens.ForwardSourceQuery()
    query.component = "thin1"
    query.min_initial_mass_msun = 0.09
    query.max_initial_mass_msun = 1.0

    # Empty selection → probability = 1
    assert generator.selection_probability(query) == 1.0

    # Wide absolute cut covering entire magnitude range → close to 1
    wide_cut = genulens.MagnitudeSelection()
    wide_cut.band = "F146mag"
    wide_cut.min_magnitude = -30.0
    wide_cut.max_magnitude = 30.0
    query.magnitude_selections = [wide_cut]
    p_wide = generator.selection_probability(query)
    assert 0.99 <= p_wide <= 1.0

    # Impossible cut (brighter than any MS star in this range)
    impossible_cut = genulens.MagnitudeSelection()
    impossible_cut.band = "F146mag"
    impossible_cut.min_magnitude = -10.0
    impossible_cut.max_magnitude = -5.0
    query.magnitude_selections = [impossible_cut]
    p_impossible = generator.selection_probability(query)
    assert p_impossible == 0.0

    # Realistic narrow cut → between wide and impossible
    narrow_cut = genulens.MagnitudeSelection()
    narrow_cut.band = "F146mag"
    narrow_cut.min_magnitude = 5.0
    narrow_cut.max_magnitude = 8.5
    query.magnitude_selections = [narrow_cut]
    p_narrow = generator.selection_probability(query)
    assert 0.0 < p_narrow < p_wide

    # Monotonicity: wider range gives higher probability
    mid_cut = genulens.MagnitudeSelection()
    mid_cut.band = "F146mag"
    mid_cut.min_magnitude = 4.0
    mid_cut.max_magnitude = 9.0
    query.magnitude_selections = [mid_cut]
    p_mid = generator.selection_probability(query)
    assert p_mid >= p_narrow


def test_source_population_prior_points():
    thin_points = genulens.SourcePopulationPrior.points_for_component(0)
    assert len(thin_points) == 12
    assert abs(sum(point.weight for point in thin_points) - 1.0) < 1e-12
    assert {point.metallicity_mh for point in thin_points} == {-0.5, -0.25, 0.0, 0.25}
    assert len({point.log_age for point in thin_points}) == 3

    thick_points = genulens.SourcePopulationPrior.points_for_component(7)
    assert len(thick_points) == 5
    assert abs(sum(point.weight for point in thick_points) - 1.0) < 1e-12
    assert -0.8 in {point.metallicity_mh for point in thick_points}


def test_selection_probability_matches_sampling():
    """selection_probability must match the IMF fraction satisfying the magnitude cut."""
    generator = genulens.ForwardSourceGenerator.load_default_roman()

    cut = genulens.MagnitudeSelection()
    cut.band = "F146mag"
    cut.min_magnitude = 5.0
    cut.max_magnitude = 8.5

    query_with_cut = genulens.ForwardSourceQuery()
    query_with_cut.component = "thin1"
    query_with_cut.min_initial_mass_msun = 0.09
    query_with_cut.max_initial_mass_msun = 1.0
    query_with_cut.magnitude_selections = [cut]
    p = generator.selection_probability(query_with_cut)
    assert 0.0 < p < 1.0

    # Sample from unrestricted IMF and measure empirical fraction
    query_no_cut = genulens.ForwardSourceQuery()
    query_no_cut.component = "thin1"
    query_no_cut.min_initial_mass_msun = 0.09
    query_no_cut.max_initial_mass_msun = 1.0
    n_total = 2000
    result = generator.sample_many(query_no_cut, n_sources=n_total, seed=42)
    arr = result.to_numpy()
    f146 = arr[:, result.columns.index("M_F146mag_S")]
    empirical_fraction = np.mean((f146 >= cut.min_magnitude) & (f146 <= cut.max_magnitude))

    # Allow 4% tolerance for Monte Carlo noise (sigma ~ sqrt(p*(1-p)/n) ≈ 1%)
    assert abs(empirical_fraction - p) < 0.04


def test_selection_probability_apparent_offset():
    """Apparent magnitude offset shifts the effective cut; closer distance = smaller offset."""
    generator = genulens.ForwardSourceGenerator.load_default_roman()

    query = genulens.ForwardSourceQuery()
    query.component = "thin1"
    query.min_initial_mass_msun = 0.09
    query.max_initial_mass_msun = 1.0

    # Apparent I=18-21 selection: at large DM+A offset, bright stars are included
    cut_near = genulens.MagnitudeSelection()
    cut_near.band = "F146mag"
    cut_near.min_magnitude = 14.0
    cut_near.max_magnitude = 22.0
    cut_near.magnitude_offset = 10.0  # small DM (nearby)
    query.magnitude_selections = [cut_near]
    p_near = generator.selection_probability(query)

    cut_far = genulens.MagnitudeSelection()
    cut_far.band = "F146mag"
    cut_far.min_magnitude = 14.0
    cut_far.max_magnitude = 22.0
    cut_far.magnitude_offset = 17.0  # large DM (far)
    query.magnitude_selections = [cut_far]
    p_far = generator.selection_probability(query)

    # Larger DM shifts faint limit fainter (more stars selected for this range)
    # or shifts beyond the grid. Either way the probabilities should differ.
    assert p_near != p_far
    assert 0.0 <= p_near <= 1.0
    assert 0.0 <= p_far <= 1.0
