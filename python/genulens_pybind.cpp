#include "genulens/options.hpp"
#include "genulens/model/forward_source.hpp"
#include "genulens/model/isochrone_grid.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/model/stellar_population_model.hpp"
#include "genulens/simulation/simulator.hpp"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(genulens, m)
{
    auto run_simulation = [](const genulens::GenulensConfig &cfg, py::object likelihood) {
        genulens::LikelihoodFunction fn;
        if (!likelihood.is_none()) {
            fn = [likelihood](const genulens::Event &event) {
                py::gil_scoped_acquire gil;
                return likelihood(event).cast<double>();
            };
        }
        py::gil_scoped_release release;
        return genulens::simulate(cfg, fn);
    };

    py::class_<genulens::GenulensConfig>(m, "Config")
        .def(py::init<>())
        .def(py::init([](double l, double b, long n_simu, unsigned long seed) {
                 genulens::GenulensConfig cfg;
                 cfg.l = l;
                 cfg.b = b;
                 cfg.n_simu = n_simu;
                 cfg.seed = seed;
                 return cfg;
             }),
             py::arg("l") = 1.0,
             py::arg("b") = -3.9,
             py::arg("n_simu") = 100000,
             py::arg("seed") = 12304357UL)
        .def_readwrite("seed", &genulens::GenulensConfig::seed)
        .def_readwrite("l", &genulens::GenulensConfig::l)
        .def_readwrite("b", &genulens::GenulensConfig::b)
        .def_readwrite("n_simu", &genulens::GenulensConfig::n_simu)
        .def_readwrite("observed_tE", &genulens::GenulensConfig::observed_tE)
        .def_readwrite("observed_tE_error", &genulens::GenulensConfig::observed_tE_error)
        .def_readwrite("model", &genulens::GenulensConfig::model)
        .def_readwrite("observation", &genulens::GenulensConfig::observation)
        .def_readwrite("source", &genulens::GenulensConfig::source)
        .def_readwrite("forward_source", &genulens::GenulensConfig::forward_source)
        .def_readwrite("sampling", &genulens::GenulensConfig::sampling)
        .def_readwrite("runtime", &genulens::GenulensConfig::runtime)
        .def_readwrite("input_dir", &genulens::GenulensConfig::input_dir);

    py::class_<genulens::SourceSelectionConfig>(m, "SourceSelectionConfig")
        .def(py::init<>())
        .def_readwrite("i_min", &genulens::SourceSelectionConfig::i_min)
        .def_readwrite("i_max", &genulens::SourceSelectionConfig::i_max)
        .def_readwrite("vi_min", &genulens::SourceSelectionConfig::vi_min)
        .def_readwrite("vi_max", &genulens::SourceSelectionConfig::vi_max)
        .def_readwrite("ai_rc", &genulens::SourceSelectionConfig::ai_rc)
        .def_readwrite("evi_rc", &genulens::SourceSelectionConfig::evi_rc)
        .def_readwrite("dm_rc", &genulens::SourceSelectionConfig::dm_rc)
        .def_readwrite("ak_rc", &genulens::SourceSelectionConfig::ak_rc)
        .def_readwrite("dust_scale_height_pc", &genulens::SourceSelectionConfig::dust_scale_height_pc);

    py::class_<genulens::ForwardSourceConfig>(m, "ForwardSourceConfig")
        .def(py::init<>())
        .def_readwrite("enabled", &genulens::ForwardSourceConfig::enabled)
        .def_readwrite("photometry", &genulens::ForwardSourceConfig::photometry)
        .def_readwrite("min_initial_mass_msun", &genulens::ForwardSourceConfig::min_initial_mass_msun)
        .def_readwrite("max_initial_mass_msun", &genulens::ForwardSourceConfig::max_initial_mass_msun);

    py::class_<genulens::SamplingConfig>(m, "SamplingConfig")
        .def(py::init<>())
        .def_readwrite("n_like_min", &genulens::SamplingConfig::n_like_min)
        .def_readwrite("v_earth_l", &genulens::SamplingConfig::v_earth_l)
        .def_readwrite("v_earth_b", &genulens::SamplingConfig::v_earth_b)
        .def_readwrite("gamma_ds", &genulens::SamplingConfig::gamma_ds)
        .def_readwrite("weight_lens_distance", &genulens::SamplingConfig::weight_lens_distance)
        .def_readwrite("weight_lens_mass", &genulens::SamplingConfig::weight_lens_mass)
        .def_readwrite("no_gamma_importance_sampling", &genulens::SamplingConfig::no_gamma_importance_sampling)
        .def_readwrite("small_gamma", &genulens::SamplingConfig::small_gamma)
        .def_readwrite("verbosity", &genulens::SamplingConfig::verbosity)
        .def_readwrite("uniform_likelihood", &genulens::SamplingConfig::uniform_likelihood)
        .def_readwrite("binary", &genulens::SamplingConfig::binary)
        .def_readwrite("remnant", &genulens::SamplingConfig::remnant)
        .def_readwrite("only_white_dwarf", &genulens::SamplingConfig::only_white_dwarf)
        .def_readwrite("calc_prior_piE", &genulens::SamplingConfig::calc_prior_piE)
        .def_readwrite("calc_prior_thetaE", &genulens::SamplingConfig::calc_prior_thetaE);

    py::class_<genulens::RuntimeConfig>(m, "RuntimeConfig")
        .def(py::init<>())
        .def_readwrite("max_distance_pc", &genulens::RuntimeConfig::max_distance_pc)
        .def_readwrite("calculate_optical_depth", &genulens::RuntimeConfig::calculate_optical_depth);

    py::class_<genulens::ObservationConfig>(m, "ObservationConfig")
        .def(py::init<>())
        .def_readwrite("tE_obs", &genulens::ObservationConfig::tE_obs)
        .def_readwrite("tE_err", &genulens::ObservationConfig::tE_err)
        .def_readwrite("tE_det", &genulens::ObservationConfig::tE_det)
        .def_readwrite("tE_min", &genulens::ObservationConfig::tE_min)
        .def_readwrite("tE_max", &genulens::ObservationConfig::tE_max)
        .def_readwrite("thetaE_obs", &genulens::ObservationConfig::thetaE_obs)
        .def_readwrite("thetaE_err", &genulens::ObservationConfig::thetaE_err)
        .def_readwrite("thetaE_det", &genulens::ObservationConfig::thetaE_det)
        .def_readwrite("thetaE_min", &genulens::ObservationConfig::thetaE_min)
        .def_readwrite("thetaE_max", &genulens::ObservationConfig::thetaE_max)
        .def_readwrite("piE_obs", &genulens::ObservationConfig::piE_obs)
        .def_readwrite("piE_err", &genulens::ObservationConfig::piE_err)
        .def_readwrite("piE_det", &genulens::ObservationConfig::piE_det)
        .def_readwrite("piE_min", &genulens::ObservationConfig::piE_min)
        .def_readwrite("piE_max", &genulens::ObservationConfig::piE_max)
        .def_readwrite("piEN_obs", &genulens::ObservationConfig::piEN_obs)
        .def_readwrite("piEN_err", &genulens::ObservationConfig::piEN_err)
        .def_readwrite("piEE_obs", &genulens::ObservationConfig::piEE_obs)
        .def_readwrite("piEE_err", &genulens::ObservationConfig::piEE_err)
        .def_readwrite("musl_obs", &genulens::ObservationConfig::musl_obs)
        .def_readwrite("musl_err", &genulens::ObservationConfig::musl_err)
        .def_readwrite("musb_obs", &genulens::ObservationConfig::musb_obs)
        .def_readwrite("musb_err", &genulens::ObservationConfig::musb_err)
        .def_readwrite("musN_obs", &genulens::ObservationConfig::musN_obs)
        .def_readwrite("musN_err", &genulens::ObservationConfig::musN_err)
        .def_readwrite("musE_obs", &genulens::ObservationConfig::musE_obs)
        .def_readwrite("musE_err", &genulens::ObservationConfig::musE_err)
        .def_readwrite("musRCG", &genulens::ObservationConfig::musRCG)
        .def_readwrite("muhelN_obs", &genulens::ObservationConfig::muhelN_obs)
        .def_readwrite("muhelN_err", &genulens::ObservationConfig::muhelN_err)
        .def_readwrite("muhelE_obs", &genulens::ObservationConfig::muhelE_obs)
        .def_readwrite("muhelE_err", &genulens::ObservationConfig::muhelE_err)
        .def_readwrite("IL_obs", &genulens::ObservationConfig::IL_obs)
        .def_readwrite("IL_err", &genulens::ObservationConfig::IL_err)
        .def_readwrite("IL_det", &genulens::ObservationConfig::IL_det)
        .def_readwrite("KL_obs", &genulens::ObservationConfig::KL_obs)
        .def_readwrite("KL_err", &genulens::ObservationConfig::KL_err)
        .def_readwrite("KL_det", &genulens::ObservationConfig::KL_det)
        .def_readwrite("u0_obs", &genulens::ObservationConfig::u0_obs);

    py::class_<genulens::model::IMFParameters>(m, "IMFParameters")
        .def(py::init<>())
        .def_readwrite("m0", &genulens::model::IMFParameters::m0)
        .def_readwrite("m1", &genulens::model::IMFParameters::m1)
        .def_readwrite("m2", &genulens::model::IMFParameters::m2)
        .def_readwrite("m3", &genulens::model::IMFParameters::m3)
        .def_readwrite("ml", &genulens::model::IMFParameters::ml)
        .def_readwrite("mu", &genulens::model::IMFParameters::mu)
        .def_readwrite("alpha0", &genulens::model::IMFParameters::alpha0)
        .def_readwrite("alpha1", &genulens::model::IMFParameters::alpha1)
        .def_readwrite("alpha2", &genulens::model::IMFParameters::alpha2)
        .def_readwrite("alpha3", &genulens::model::IMFParameters::alpha3)
        .def_readwrite("alpha4", &genulens::model::IMFParameters::alpha4);

    py::class_<genulens::model::DensityParameters>(m, "DensityParameters")
        .def(py::init<>())
        .def_readwrite("disk", &genulens::model::DensityParameters::disk)
        .def_readwrite("rho_t0", &genulens::model::DensityParameters::rho_t0)
        .def_readwrite("h_disk", &genulens::model::DensityParameters::h_disk)
        .def_readwrite("add_x", &genulens::model::DensityParameters::add_x)
        .def_readwrite("model", &genulens::model::DensityParameters::model)
        .def_readwrite("r0", &genulens::model::DensityParameters::r0)
        .def_readwrite("theta_d", &genulens::model::DensityParameters::theta_d)
        .def_readwrite("frho0b", &genulens::model::DensityParameters::frho0b)
        .def_readwrite("rc", &genulens::model::DensityParameters::rc)
        .def_readwrite("zb_c", &genulens::model::DensityParameters::zb_c)
        .def_readwrite("x0", &genulens::model::DensityParameters::x0)
        .def_readwrite("y0", &genulens::model::DensityParameters::y0)
        .def_readwrite("z0", &genulens::model::DensityParameters::z0)
        .def_readwrite("c1", &genulens::model::DensityParameters::c1)
        .def_readwrite("c2", &genulens::model::DensityParameters::c2)
        .def_readwrite("c3", &genulens::model::DensityParameters::c3)
        .def_readwrite("x0_x", &genulens::model::DensityParameters::x0_x)
        .def_readwrite("y0_x", &genulens::model::DensityParameters::y0_x)
        .def_readwrite("z0_x", &genulens::model::DensityParameters::z0_x)
        .def_readwrite("c1_x", &genulens::model::DensityParameters::c1_x)
        .def_readwrite("c2_x", &genulens::model::DensityParameters::c2_x)
        .def_readwrite("b_zx", &genulens::model::DensityParameters::b_zx)
        .def_readwrite("f_x", &genulens::model::DensityParameters::f_x)
        .def_readwrite("rc_x", &genulens::model::DensityParameters::rc_x)
        .def_readwrite("b_zy", &genulens::model::DensityParameters::b_zy)
        .def_readwrite("stellar_halo", &genulens::model::DensityParameters::stellar_halo)
        .def_readwrite("rho0_sh_ms", &genulens::model::DensityParameters::rho0_sh_ms);

    py::class_<genulens::model::KinematicsParameters>(m, "KinematicsParameters")
        .def(py::init<>())
        .def_readwrite("omega_p", &genulens::model::KinematicsParameters::omega_p)
        .def_readwrite("model_vb", &genulens::model::KinematicsParameters::model_vb)
        .def_readwrite("x0_vb", &genulens::model::KinematicsParameters::x0_vb)
        .def_readwrite("y0_vb", &genulens::model::KinematicsParameters::y0_vb)
        .def_readwrite("z0_vb", &genulens::model::KinematicsParameters::z0_vb)
        .def_readwrite("c1_vb", &genulens::model::KinematicsParameters::c1_vb)
        .def_readwrite("c2_vb", &genulens::model::KinematicsParameters::c2_vb)
        .def_readwrite("c3_vb", &genulens::model::KinematicsParameters::c3_vb)
        .def_readwrite("sigx_vb", &genulens::model::KinematicsParameters::sigx_vb)
        .def_readwrite("sigy_vb", &genulens::model::KinematicsParameters::sigy_vb)
        .def_readwrite("sigz_vb", &genulens::model::KinematicsParameters::sigz_vb)
        .def_readwrite("sigx_vb0", &genulens::model::KinematicsParameters::sigx_vb0)
        .def_readwrite("sigy_vb0", &genulens::model::KinematicsParameters::sigy_vb0)
        .def_readwrite("sigz_vb0", &genulens::model::KinematicsParameters::sigz_vb0)
        .def_readwrite("vx_stream", &genulens::model::KinematicsParameters::vx_stream)
        .def_readwrite("y0_stream", &genulens::model::KinematicsParameters::y0_stream)
        .def_readwrite("model_vbz", &genulens::model::KinematicsParameters::model_vbz)
        .def_readwrite("x0_vbz", &genulens::model::KinematicsParameters::x0_vbz)
        .def_readwrite("y0_vbz", &genulens::model::KinematicsParameters::y0_vbz)
        .def_readwrite("z0_vbz", &genulens::model::KinematicsParameters::z0_vbz)
        .def_readwrite("c1_vbz", &genulens::model::KinematicsParameters::c1_vbz)
        .def_readwrite("c2_vbz", &genulens::model::KinematicsParameters::c2_vbz)
        .def_readwrite("c3_vbz", &genulens::model::KinematicsParameters::c3_vbz)
        .def_readwrite("hsig_ut", &genulens::model::KinematicsParameters::hsig_ut)
        .def_readwrite("hsig_wt", &genulens::model::KinematicsParameters::hsig_wt)
        .def_readwrite("hsig_uT", &genulens::model::KinematicsParameters::hsig_uT)
        .def_readwrite("hsig_wT", &genulens::model::KinematicsParameters::hsig_wT)
        .def_readwrite("beta_u", &genulens::model::KinematicsParameters::beta_u)
        .def_readwrite("beta_w", &genulens::model::KinematicsParameters::beta_w)
        .def_readwrite("sig_u10d", &genulens::model::KinematicsParameters::sig_u10d)
        .def_readwrite("sig_w10d", &genulens::model::KinematicsParameters::sig_w10d)
        .def_readwrite("sig_u0td", &genulens::model::KinematicsParameters::sig_u0td)
        .def_readwrite("sig_w0td", &genulens::model::KinematicsParameters::sig_w0td)
        .def_readwrite("sig_u_sh", &genulens::model::KinematicsParameters::sig_u_sh)
        .def_readwrite("sig_v_sh", &genulens::model::KinematicsParameters::sig_v_sh)
        .def_readwrite("sig_w_sh", &genulens::model::KinematicsParameters::sig_w_sh);

    py::class_<genulens::model::NsdParameters>(m, "NsdParameters")
        .def(py::init<>())
        .def_readwrite("enabled", &genulens::model::NsdParameters::enabled)
        .def_readwrite("x0", &genulens::model::NsdParameters::x0)
        .def_readwrite("y0", &genulens::model::NsdParameters::y0)
        .def_readwrite("z0", &genulens::model::NsdParameters::z0)
        .def_readwrite("mass", &genulens::model::NsdParameters::mass);

    py::class_<genulens::model::BhKickParameters>(m, "BhKickParameters")
        .def(py::init<>())
        .def_readwrite("mix_disk_kick", &genulens::model::BhKickParameters::mix_disk_kick)
        .def_readwrite("kick_bh", &genulens::model::BhKickParameters::kick_bh)
        .def_readwrite("kick_ns", &genulens::model::BhKickParameters::kick_ns)
        .def_readwrite("disk_scale_height", &genulens::model::BhKickParameters::disk_scale_height)
        .def_readwrite("bar_scale_height", &genulens::model::BhKickParameters::bar_scale_height)
        .def_readwrite("fix_disk_scale_length", &genulens::model::BhKickParameters::fix_disk_scale_length)
        .def_readwrite("disk_scale_length", &genulens::model::BhKickParameters::disk_scale_length)
        .def_readwrite("beta", &genulens::model::BhKickParameters::beta)
        .def_readwrite("use_sigma_correction", &genulens::model::BhKickParameters::use_sigma_correction);

    py::class_<genulens::model::ModelParameters>(m, "ModelParameters")
        .def(py::init<>())
        .def_readwrite("imf", &genulens::model::ModelParameters::imf)
        .def_readwrite("density", &genulens::model::ModelParameters::density)
        .def_readwrite("kinematics", &genulens::model::ModelParameters::kinematics)
        .def_readwrite("nsd", &genulens::model::ModelParameters::nsd)
        .def_readwrite("bh_kick", &genulens::model::ModelParameters::bh_kick);

    py::class_<genulens::model::IsochroneQuery>(m, "IsochroneQuery")
        .def(py::init<>())
        .def_readwrite("log_age", &genulens::model::IsochroneQuery::log_age)
        .def_readwrite("metallicity_mh", &genulens::model::IsochroneQuery::metallicity_mh)
        .def_readwrite("initial_mass_msun", &genulens::model::IsochroneQuery::initial_mass_msun)
        .def_readwrite("component", &genulens::model::IsochroneQuery::component);

    py::class_<genulens::model::StellarProperties>(m, "StellarProperties")
        .def_readonly("component", &genulens::model::StellarProperties::component)
        .def_readonly("component_index", &genulens::model::StellarProperties::component_index)
        .def_readonly("family", &genulens::model::StellarProperties::family)
        .def_readonly("log_age", &genulens::model::StellarProperties::log_age)
        .def_readonly("metallicity_mh", &genulens::model::StellarProperties::metallicity_mh)
        .def_readonly("zini", &genulens::model::StellarProperties::zini)
        .def_readonly("initial_mass_msun", &genulens::model::StellarProperties::initial_mass_msun)
        .def_readonly("current_mass_msun", &genulens::model::StellarProperties::current_mass_msun)
        .def_readonly("radius_rsun", &genulens::model::StellarProperties::radius_rsun)
        .def_readonly("teff_k", &genulens::model::StellarProperties::teff_k)
        .def_readonly("logg", &genulens::model::StellarProperties::logg)
        .def_readonly("absolute_magnitudes", &genulens::model::StellarProperties::absolute_magnitudes);

    py::class_<genulens::model::IsochroneGrid>(m, "IsochroneGrid")
        .def_static("load", &genulens::model::IsochroneGrid::load, py::arg("path"))
        .def_static("load_default_roman", &genulens::model::IsochroneGrid::load_default_roman)
        .def_static("load_default_prime", &genulens::model::IsochroneGrid::load_default_prime)
        .def_property_readonly("bands", &genulens::model::IsochroneGrid::bands)
        .def_property_readonly("row_count", &genulens::model::IsochroneGrid::row_count)
        .def_property_readonly("sequence_count", &genulens::model::IsochroneGrid::sequence_count)
        .def("lookup", &genulens::model::IsochroneGrid::lookup, py::arg("query"));

    py::class_<genulens::model::StellarPopulationQuery>(m, "StellarPopulationQuery")
        .def(py::init<>())
        .def_readwrite("component", &genulens::model::StellarPopulationQuery::component)
        .def_readwrite("component_index", &genulens::model::StellarPopulationQuery::component_index)
        .def_readwrite("initial_mass_msun", &genulens::model::StellarPopulationQuery::initial_mass_msun)
        .def_readwrite("log_age", &genulens::model::StellarPopulationQuery::log_age)
        .def_readwrite("metallicity_mh", &genulens::model::StellarPopulationQuery::metallicity_mh)
        .def_readwrite("use_default_log_age", &genulens::model::StellarPopulationQuery::use_default_log_age)
        .def_readwrite("use_default_metallicity", &genulens::model::StellarPopulationQuery::use_default_metallicity);

    py::class_<genulens::model::StellarPopulationModel>(m, "StellarPopulationModel")
        .def_static("load_default_roman", &genulens::model::StellarPopulationModel::load_default_roman)
        .def_static("load_default_prime", &genulens::model::StellarPopulationModel::load_default_prime)
        .def("lookup", &genulens::model::StellarPopulationModel::lookup, py::arg("query"))
        .def_static("component_name", &genulens::model::StellarPopulationModel::component_name, py::arg("component_index"))
        .def_static("component_index", &genulens::model::StellarPopulationModel::component_index, py::arg("component"))
        .def_static("default_log_age", &genulens::model::StellarPopulationModel::default_log_age,
                    py::arg("component"))
        .def_static("default_metallicity_mh", &genulens::model::StellarPopulationModel::default_metallicity_mh,
                    py::arg("component"));

    py::class_<genulens::model::ForwardSourceQuery>(m, "ForwardSourceQuery")
        .def(py::init<>())
        .def_readwrite("component", &genulens::model::ForwardSourceQuery::component)
        .def_readwrite("component_index", &genulens::model::ForwardSourceQuery::component_index)
        .def_readwrite("distance_pc", &genulens::model::ForwardSourceQuery::distance_pc)
        .def_readwrite("min_initial_mass_msun", &genulens::model::ForwardSourceQuery::min_initial_mass_msun)
        .def_readwrite("max_initial_mass_msun", &genulens::model::ForwardSourceQuery::max_initial_mass_msun)
        .def_readwrite("log_age", &genulens::model::ForwardSourceQuery::log_age)
        .def_readwrite("metallicity_mh", &genulens::model::ForwardSourceQuery::metallicity_mh)
        .def_readwrite("use_default_log_age", &genulens::model::ForwardSourceQuery::use_default_log_age)
        .def_readwrite("use_default_metallicity", &genulens::model::ForwardSourceQuery::use_default_metallicity);

    py::class_<genulens::model::ForwardSource>(m, "ForwardSource")
        .def_readonly("stellar", &genulens::model::ForwardSource::stellar)
        .def_readonly("distance_pc", &genulens::model::ForwardSource::distance_pc)
        .def_readonly("angular_radius_microarcsec", &genulens::model::ForwardSource::angular_radius_microarcsec);

    py::class_<genulens::model::ForwardSourceResult>(m, "ForwardSourceResult")
        .def_property_readonly("columns", &genulens::model::ForwardSourceResult::columns)
        .def_property_readonly("bands", [](const genulens::model::ForwardSourceResult &result) {
            return result.bands;
        })
        .def_property_readonly("row_count", &genulens::model::ForwardSourceResult::row_count)
        .def("to_numpy", [](const genulens::model::ForwardSourceResult &result) {
            auto rows = result.flattened_rows();
            py::array_t<double> array({result.row_count(), result.column_count()});
            auto mutable_array = array.mutable_unchecked<2>();
            for (py::ssize_t r = 0; r < static_cast<py::ssize_t>(result.row_count()); ++r) {
                for (py::ssize_t c = 0; c < static_cast<py::ssize_t>(result.column_count()); ++c) {
                    mutable_array(r, c) = rows[static_cast<std::size_t>(r * result.column_count() + c)];
                }
            }
            return array;
        });

    py::class_<genulens::model::ForwardSourceGenerator>(m, "ForwardSourceGenerator")
        .def_static("load_default_roman", &genulens::model::ForwardSourceGenerator::load_default_roman,
                    py::arg("imf_parameters") = genulens::model::default_model_parameters().imf)
        .def_static("load_default_prime", &genulens::model::ForwardSourceGenerator::load_default_prime,
                    py::arg("imf_parameters") = genulens::model::default_model_parameters().imf)
        .def("sample", [](const genulens::model::ForwardSourceGenerator &generator,
                          const genulens::model::ForwardSourceQuery &query,
                          unsigned long seed) {
            genulens::RandomEngine rng(seed);
            return generator.sample(query, rng);
        }, py::arg("query"), py::arg("seed") = 12304357UL)
        .def("sample_many", [](const genulens::model::ForwardSourceGenerator &generator,
                               const genulens::model::ForwardSourceQuery &query,
                               std::size_t n_sources,
                               unsigned long seed) {
            genulens::RandomEngine rng(seed);
            return generator.sample_many(query, n_sources, rng);
        }, py::arg("query"), py::arg("n_sources"), py::arg("seed") = 12304357UL);
    m.def("angular_radius_microarcsec", &genulens::model::angular_radius_microarcsec,
          py::arg("radius_rsun"), py::arg("distance_pc"));

    py::class_<genulens::Event>(m, "Event")
        .def_readonly("wtj", &genulens::Event::weight)
        .def_readonly("t_E", &genulens::Event::tE)
        .def_readonly("theta_E", &genulens::Event::thetaE)
        .def_readonly("pi_E", &genulens::Event::piE)
        .def_readonly("pi_EN", &genulens::Event::piEN)
        .def_readonly("pi_EE", &genulens::Event::piEE)
        .def_readonly("D_L", &genulens::Event::lens_distance_pc)
        .def_readonly("D_S", &genulens::Event::source_distance_pc)
        .def_readonly("M_L", &genulens::Event::lens_mass_msun)
        .def_readonly("mu_rel", &genulens::Event::mu_rel_masyr)
        .def_readonly("mu_rel_N", &genulens::Event::mu_rel_N_masyr)
        .def_readonly("mu_rel_E", &genulens::Event::mu_rel_E_masyr)
        .def_readonly("mu_Sl", &genulens::Event::source_mu_l_masyr)
        .def_readonly("mu_Sb", &genulens::Event::source_mu_b_masyr)
        .def_readonly("I_L", &genulens::Event::lens_i_mag)
        .def_readonly("K_L", &genulens::Event::lens_k_mag)
        .def_readonly("iL", &genulens::Event::lens_component)
        .def_readonly("iS", &genulens::Event::source_component)
        .def_readonly("fREM", &genulens::Event::remnant_flag)
        .def_readonly("tau_s", &genulens::Event::source_age_gyr)
        .def_readonly("tau_l", &genulens::Event::lens_age_gyr)
        .def_readonly("logage_S", &genulens::Event::source_log_age)
        .def_readonly("MH_S", &genulens::Event::source_metallicity_mh)
        .def_readonly("M_S_ini", &genulens::Event::source_initial_mass_msun)
        .def_readonly("M_S", &genulens::Event::source_current_mass_msun)
        .def_readonly("R_S", &genulens::Event::source_radius_rsun)
        .def_readonly("teff_S", &genulens::Event::source_teff_k)
        .def_readonly("logg_S", &genulens::Event::source_logg)
        .def_readonly("theta_S", &genulens::Event::source_angular_radius_microarcsec)
        .def_readonly("Mabs_S", &genulens::Event::source_absolute_magnitudes);

    py::class_<genulens::SimulationResult>(m, "SimulationResult")
        .def_property_readonly("columns", &genulens::SimulationResult::columns)
        .def("to_numpy", [](const genulens::SimulationResult &result) {
            auto rows = result.flattened_rows();
            py::array_t<double> array({result.row_count(), result.column_count()});
            auto mutable_array = array.mutable_unchecked<2>();
            for (py::ssize_t r = 0; r < static_cast<py::ssize_t>(result.row_count()); ++r) {
                for (py::ssize_t c = 0; c < static_cast<py::ssize_t>(result.column_count()); ++c) {
                    mutable_array(r, c) = rows[static_cast<std::size_t>(r * result.column_count() + c)];
                }
            }
            return array;
        });

    m.def("simulate", run_simulation, py::arg("cfg"), py::arg("likelihood") = py::none());

    m.def("ruc", [run_simulation](double l, double b, long n_simu, unsigned long seed, py::object likelihood) {
        genulens::GenulensConfig cfg;
        cfg.l = l;
        cfg.b = b;
        cfg.n_simu = n_simu;
        cfg.seed = seed;
        return run_simulation(cfg, likelihood);
    }, py::arg("l") = 1.0, py::arg("b") = -3.9, py::arg("n_simu") = 100000,
       py::arg("seed") = 12304357UL, py::arg("likelihood") = py::none());
}
