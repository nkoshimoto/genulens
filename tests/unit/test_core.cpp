#include "genulens/io/input_data.hpp"
#include "genulens/math/quadrature.hpp"
#include "genulens/model/coordinates.hpp"
#include "genulens/model/kinematics.hpp"
#include "genulens/model/mass_function.hpp"
#include "genulens/model/parameters.hpp"
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

    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(std::filesystem::temp_directory_path());
    try {
        require(genulens::resolve_input_file("Minidie.dat").filename() == "Minidie.dat",
                "input file resolution from another cwd failed");
        require(genulens::open_input_file("input_files/Minidie.dat", "r") != nullptr,
                "input fopen compatibility from another cwd failed");
    } catch (...) {
        std::filesystem::current_path(cwd);
        throw;
    }
    std::filesystem::current_path(cwd);

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
