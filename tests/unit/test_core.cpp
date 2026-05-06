#include "genulens/io/input_data.hpp"
#include "genulens/math/quadrature.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/forward_source.hpp"
#include "genulens/model/isochrone_grid.hpp"
#include "genulens/model/kinematics.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"
#include "genulens/model/stellar_population_model.hpp"
#include "genulens/rng.hpp"
#include "genulens/simulation/likelihood.hpp"
#include "genulens/simulation/observation_likelihood.hpp"
#include "genulens/simulation/simulator.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main()
{
    genulens::RandomEngine rng1(1234);
    genulens::RandomEngine rng2(1234);
    require(rng1.uniform() == rng2.uniform(), "RandomEngine seed reproducibility failed");

    const auto xyz = genulens::model::galactic_to_cartesian(0.0, {1.0, -3.9});
    require(std::abs(xyz.x - 8160.0) < 1e-9, "coordinate conversion at Sun failed");
    const auto pa = genulens::model::CoordinateTransformer::position_angle(1.0, -3.9);
    require(std::isfinite(pa.degrees), "position angle calculation failed");

    double locations[64] = {};
    double weights[64] = {};
    const int nmin = genulens::math::NewtonCotes::coefficients(4, locations, weights);
    require(nmin == 7, "Newton-Cotes order normalization failed");
    require(std::abs(locations[1] - 0.25) < 1e-12, "Newton-Cotes locations changed");
    require(std::abs(weights[0] - 70.0 / 360.0) < 1e-12, "Newton-Cotes weights changed");

    const auto grid = genulens::model::BrokenPowerLawIMF(genulens::model::default_model_parameters().imf)
                          .build_grid(100, 0.001, 120.0);
    require(grid.log_mass.size() == 101, "IMF grid size changed");
    require(std::abs(grid.cumulative_number_norm.front()) < 1e-15, "IMF cumulative start changed");
    require(std::abs(grid.cumulative_number_norm.back() - 1.0) < 1e-12, "IMF cumulative normalization changed");
    genulens::RandomEngine remnant_rng(10);
    const auto wd = genulens::model::RemnantMassModel(9.0).evolve(1.0, true, remnant_rng);
    require(std::abs(wd.mass_msun - 0.503) < 1e-12, "WD remnant mass relation changed");
    genulens::RandomEngine binary_rng(10);
    const auto separation = genulens::model::BinaryLensSampler().sample_projected_separation(0.5, 0.2, 1, binary_rng);
    require(std::isfinite(separation.log_separation_au), "binary log separation failed");
    require(separation.projected_separation_au > 0.0, "binary projected separation failed");

    genulens::Event event;
    event.tE = 10.0;
    genulens::GaussianLikelihood like(10.0, 1.0);
    require(std::abs(like(event) - 1.0) < 1e-12, "Gaussian likelihood peak failed");

    genulens::RandomEngine like_rng(42);
    genulens::ObservationConstraint constraint;
    constraint.observed = 10.0;
    constraint.error = 1.0;
    constraint.uniform = true;
    require(genulens::ObservationLikelihood(constraint).accept_weight(10.5, like_rng) == 1.0,
            "uniform observation likelihood accept failed");

    genulens::InputDataRepository repo;
    require(repo.resolve("Minidie.dat").filename() == "Minidie.dat", "input file resolution failed");
    require(repo.resolve("input_files/Minidie.dat").filename() == "Minidie.dat", "prefixed input file resolution failed");
    require(repo.resolve("source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat").filename() ==
                "all_roman_parsec.dat",
            "nested input file resolution failed");

    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(std::filesystem::temp_directory_path());
    try {
        require(genulens::resolve_input_file("Minidie.dat").filename() == "Minidie.dat",
                "input file resolution from another cwd failed");
        require(genulens::resolve_input_file("input_files/source_photometry/parsec_cmd/metallicity_grid/normalized/all_roman_parsec.dat")
                    .filename() == "all_roman_parsec.dat",
                "prefixed nested input file resolution from another cwd failed");
        require(genulens::open_input_file("input_files/Minidie.dat", "r") != nullptr,
                "input fopen compatibility from another cwd failed");
    } catch (...) {
        std::filesystem::current_path(cwd);
        throw;
    }
    std::filesystem::current_path(cwd);

    const auto roman_isochrones = genulens::model::IsochroneGrid::load_default_roman();
    require(roman_isochrones.row_count() == 14723, "Roman isochrone row count changed");
    require(roman_isochrones.sequence_count() == 41, "Roman isochrone sequence count changed");
    require(roman_isochrones.bands().size() == 6, "Roman isochrone band count changed");

    genulens::model::IsochroneQuery iso_query;
    iso_query.component = "thin1";
    iso_query.log_age = 8.0;
    iso_query.metallicity_mh = -0.5;
    iso_query.initial_mass_msun = 0.1000000015;
    const auto roman_star = roman_isochrones.lookup(iso_query);
    require(roman_star.component == "thin1", "isochrone component lookup failed");
    require(std::abs(roman_star.metallicity_mh + 0.5) < 1e-12, "isochrone metallicity lookup failed");
    require(std::abs(roman_star.initial_mass_msun - 0.1000000015) < 1e-9, "isochrone mass lookup failed");
    require(std::abs(roman_star.teff_k - 2886.02441) < 1e-5, "isochrone Teff lookup failed");
    require(roman_star.absolute_magnitudes.count("F146mag") == 1, "Roman isochrone F146 band missing");
    require(std::abs(roman_star.absolute_magnitudes.at("F146mag") - 9.251) < 1e-12,
            "Roman isochrone F146 lookup failed");

    iso_query.metallicity_mh = -0.45;
    const auto nearest_mh_star = roman_isochrones.lookup(iso_query);
    require(std::abs(nearest_mh_star.metallicity_mh + 0.5) < 1e-12,
            "isochrone nearest metallicity lookup failed");

    const auto prime_isochrones = genulens::model::IsochroneGrid::load_default_prime();
    const auto prime_star = prime_isochrones.lookup(iso_query);
    require(prime_star.absolute_magnitudes.count("Vmag") == 1, "prime isochrone V band missing");
    require(prime_star.absolute_magnitudes.count("Imag") == 1, "prime isochrone I band missing");

    const auto population = genulens::model::StellarPopulationModel::load_default_roman();
    genulens::model::StellarPopulationQuery pop_query;
    pop_query.component_index = 7;
    pop_query.initial_mass_msun = 0.1000000015;
    const auto thick_star = population.lookup(pop_query);
    require(thick_star.component == "thick", "population component-index lookup failed");
    require(std::abs(thick_star.log_age - 10.07918) < 1e-12, "population default thick age failed");
    require(std::abs(thick_star.metallicity_mh + 0.8) < 1e-12, "population default thick metallicity failed");

    pop_query.component = "thin1";
    pop_query.component_index = -1;
    const auto thin_star = population.lookup(pop_query);
    require(thin_star.component == "thin1", "population component-name lookup failed");
    require(std::abs(thin_star.log_age - 8.0) < 1e-12, "population default thin age failed");
    require(std::abs(thin_star.metallicity_mh) < 1e-12, "population default thin metallicity failed");

    pop_query.use_default_metallicity = false;
    pop_query.metallicity_mh = -0.5;
    const auto metal_poor_thin = population.lookup(pop_query);
    require(std::abs(metal_poor_thin.metallicity_mh + 0.5) < 1e-12,
            "population explicit metallicity failed");

    genulens::RandomEngine source_rng(123);
    const auto source_generator = genulens::model::ForwardSourceGenerator::load_default_roman();
    genulens::model::ForwardSourceQuery source_query;
    source_query.component = "thin1";
    source_query.distance_pc = 8000.0;
    source_query.min_initial_mass_msun = 0.1;
    source_query.max_initial_mass_msun = 0.11;
    const auto source_star = source_generator.sample(source_query, source_rng);
    require(source_star.stellar.component == "thin1", "forward source component failed");
    require(source_star.stellar.initial_mass_msun >= 0.1 &&
                source_star.stellar.initial_mass_msun <= 0.11,
            "forward source mass range failed");
    require(source_star.angular_radius_microarcsec > 0.0, "forward source angular radius failed");
    require(source_star.apparent_magnitudes.count("F146mag") == 1, "forward source apparent F146 missing");
    require(std::abs(source_star.apparent_magnitudes.at("F146mag") -
                     source_star.stellar.absolute_magnitudes.at("F146mag") -
                     genulens::model::distance_modulus(8000.0)) < 1e-12,
            "forward source apparent magnitude failed");

    genulens::GenulensConfig cfg;
    cfg.n_simu = 5;
    const auto result = genulens::simulate(cfg, [](const genulens::Event &e) {
        return e.tE > 0.0 ? 1.0 : 0.0;
    });
    require(result.events.size() == 5, "simulator result size failed");
    require(result.column_count() == 10, "result column count failed");

    std::cout << "core unit tests passed\n";
    return 0;
}
