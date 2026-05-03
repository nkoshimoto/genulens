cmake_minimum_required(VERSION 3.16)
project(genulens LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)

option(GENULENS_BUILD_PYTHON "Build pybind11 Python extension" ON)

find_package(GSL QUIET)
if(NOT GSL_FOUND)
    find_path(GSL_INCLUDE_DIR gsl/gsl_rng.h
        HINTS
            "$ENV{GSL_ROOT}/include"
            "/rogue1_8/nunota/local/gsl/include"
            "/rogue1_8/nunota/local/gulls_env/include"
            "/usr/include"
            "/usr/local/include")
    find_library(GSL_LIBRARY gsl
        HINTS
            "$ENV{GSL_ROOT}/lib"
            "/rogue1_8/nunota/local/gsl/lib"
            "/rogue1_8/nunota/local/gulls_env/lib"
            "/usr/lib64"
            "/usr/lib"
            "/usr/local/lib")
    find_library(GSLCBLAS_LIBRARY gslcblas
        HINTS
            "$ENV{GSL_ROOT}/lib"
            "/rogue1_8/nunota/local/gsl/lib"
            "/rogue1_8/nunota/local/gulls_env/lib"
            "/usr/lib64"
            "/usr/lib"
            "/usr/local/lib")
    if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY OR NOT GSLCBLAS_LIBRARY)
        message(FATAL_ERROR "GSL was not found. Set GSL_ROOT or install gsl-config/pkg-config metadata.")
    endif()
    add_library(GSL::gsl UNKNOWN IMPORTED)
    set_target_properties(GSL::gsl PROPERTIES
        IMPORTED_LOCATION "${GSL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}")
    add_library(GSL::gslcblas UNKNOWN IMPORTED)
    set_target_properties(GSL::gslcblas PROPERTIES IMPORTED_LOCATION "${GSLCBLAS_LIBRARY}")
endif()

add_library(genulens_core
    src/genulens/types.cpp
    src/genulens/rng.cpp
    src/genulens/cli/option.cpp
    src/genulens/model/parameters.cpp
    src/genulens/math/interpolation.cpp
    src/genulens/math/quadrature.cpp
    src/genulens/model/coordinates.cpp
    src/genulens/model/galactic_model.cpp
    src/genulens/model/density.cpp
    src/genulens/model/kinematics.cpp
    src/genulens/model/mass_function.cpp
    src/genulens/model/extinction.cpp
    src/genulens/simulation/initialize.cpp
    src/genulens/simulation/sampler.cpp
    src/genulens/simulation/observation_runtime.cpp
    src/genulens/simulation/stellar_population_runtime.cpp
    src/genulens/simulation/kinematics_runtime.cpp
    src/genulens/simulation/density_runtime.cpp
    src/genulens/simulation/math_runtime.cpp
    src/genulens/simulation/cli.cpp
    src/genulens/simulation/backend.cpp
    src/genulens/simulation/simulator.cpp
    src/genulens/simulation/likelihood.cpp
    src/genulens/simulation/observation_likelihood.cpp
    src/genulens/io/input_data.cpp)
target_include_directories(genulens_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(genulens_core PUBLIC GSL::gsl GSL::gslcblas m)
target_compile_definitions(genulens_core PRIVATE GENULENS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(genulens genulens.cpp)
target_link_libraries(genulens PRIVATE genulens_core)

add_executable(calc_rho_profile
    src/genulens/tools/pre_gapmoe/calc_rho_profile.cpp
    src/genulens/tools/pre_gapmoe/galactic_model.cpp
    src/genulens/tools/pre_gapmoe/galactic_kinematics.cpp
    src/genulens/tools/pre_gapmoe/option.cpp)
target_include_directories(calc_rho_profile PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/src/genulens/tools/pre_gapmoe" ${GSL_INCLUDE_DIR})
target_link_libraries(calc_rho_profile PRIVATE genulens_core GSL::gsl GSL::gslcblas m)
set_target_properties(calc_rho_profile PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pre_gapmoe")

add_executable(calc_mass_dist
    src/genulens/tools/pre_gapmoe/calc_mass_dist.cpp
    src/genulens/tools/pre_gapmoe/galactic_model.cpp
    src/genulens/tools/pre_gapmoe/galactic_kinematics.cpp
    src/genulens/tools/pre_gapmoe/option.cpp)
target_include_directories(calc_mass_dist PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/src/genulens/tools/pre_gapmoe" ${GSL_INCLUDE_DIR})
target_link_libraries(calc_mass_dist PRIVATE genulens_core GSL::gsl GSL::gslcblas m)
set_target_properties(calc_mass_dist PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pre_gapmoe")

add_executable(calc_murel_dist
    src/genulens/tools/pre_gapmoe/calc_murel_dist.cpp
    src/genulens/tools/pre_gapmoe/galactic_model.cpp
    src/genulens/tools/pre_gapmoe/galactic_kinematics.cpp
    src/genulens/tools/pre_gapmoe/option.cpp)
target_include_directories(calc_murel_dist PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/src/genulens/tools/pre_gapmoe" ${GSL_INCLUDE_DIR})
target_link_libraries(calc_murel_dist PRIVATE genulens_core GSL::gsl GSL::gslcblas m)
set_target_properties(calc_murel_dist PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pre_gapmoe")

add_executable(test_core tests/unit/test_core.cpp)
target_link_libraries(test_core PRIVATE genulens_core)
add_test(NAME unit_core COMMAND test_core)
set_tests_properties(unit_core PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME smoke_cli COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/smoke_cli.sh")
add_test(NAME regression_cli COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/regression_cli.sh")

find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)
find_package(pybind11 CONFIG QUIET)
if(GENULENS_BUILD_PYTHON AND Python3_FOUND AND pybind11_FOUND)
    pybind11_add_module(genulens_python python/genulens_pybind.cpp)
    target_link_libraries(genulens_python PRIVATE genulens_core)
    set_target_properties(genulens_python PROPERTIES OUTPUT_NAME genulens)
    add_test(NAME python_binding
        COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/test_python.py")
    set_tests_properties(python_binding PROPERTIES
        ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:genulens_python>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
else()
    message(STATUS "Python extension disabled or pybind11/Python development files unavailable")
endif()

install(TARGETS genulens genulens_core)
install(DIRECTORY input_files DESTINATION share/genulens)
