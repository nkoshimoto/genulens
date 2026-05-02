#include "genulens/options.hpp"
#include "genulens/simulation/simulator.hpp"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(genulens, m)
{
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
        .def_readwrite("input_dir", &genulens::GenulensConfig::input_dir);

    py::class_<genulens::Event>(m, "Event")
        .def_readonly("weight", &genulens::Event::weight)
        .def_readonly("tE", &genulens::Event::tE)
        .def_readonly("thetaE", &genulens::Event::thetaE)
        .def_readonly("piE", &genulens::Event::piE)
        .def_readonly("lens_distance_pc", &genulens::Event::lens_distance_pc)
        .def_readonly("source_distance_pc", &genulens::Event::source_distance_pc)
        .def_readonly("lens_mass_msun", &genulens::Event::lens_mass_msun)
        .def_readonly("mu_rel_masyr", &genulens::Event::mu_rel_masyr)
        .def_readonly("lens_component", &genulens::Event::lens_component)
        .def_readonly("source_component", &genulens::Event::source_component);

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

    m.def("simulate", [](const genulens::GenulensConfig &cfg, py::object likelihood) {
        genulens::LikelihoodFunction fn;
        if (!likelihood.is_none()) {
            fn = [likelihood](const genulens::Event &event) {
                py::gil_scoped_acquire gil;
                return likelihood(event).cast<double>();
            };
        }
        py::gil_scoped_release release;
        return genulens::simulate(cfg, fn);
    }, py::arg("cfg"), py::arg("likelihood") = py::none());
}

