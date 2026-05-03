#include "genulens/simulation/scientific_backend.hpp"

#include "genulens/simulation/scientific_cli.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace genulens {

namespace {

class StdoutCapture {
public:
    StdoutCapture()
    {
        char pattern[] = "/tmp/genulens-capture-XXXXXX";
        fd_ = mkstemp(pattern);
        if (fd_ < 0) {
            throw std::runtime_error(std::string("mkstemp failed: ") + std::strerror(errno));
        }
        path_ = pattern;
        saved_stdout_ = dup(STDOUT_FILENO);
        if (saved_stdout_ < 0) {
            cleanup();
            throw std::runtime_error(std::string("dup stdout failed: ") + std::strerror(errno));
        }
        std::fflush(stdout);
        if (dup2(fd_, STDOUT_FILENO) < 0) {
            cleanup();
            throw std::runtime_error(std::string("redirect stdout failed: ") + std::strerror(errno));
        }
    }

    ~StdoutCapture()
    {
        cleanup();
    }

    StdoutCapture(const StdoutCapture &) = delete;
    StdoutCapture &operator=(const StdoutCapture &) = delete;

    std::string finish()
    {
        if (finished_) return output_;
        std::fflush(stdout);
        if (saved_stdout_ >= 0) {
            dup2(saved_stdout_, STDOUT_FILENO);
            close(saved_stdout_);
            saved_stdout_ = -1;
        }
        if (fd_ >= 0) {
            close(fd_);
            fd_ = -1;
        }
        std::ifstream in(path_);
        std::ostringstream buffer;
        buffer << in.rdbuf();
        output_ = buffer.str();
        unlink(path_.c_str());
        finished_ = true;
        return output_;
    }

private:
    void cleanup()
    {
        if (saved_stdout_ >= 0) {
            std::fflush(stdout);
            dup2(saved_stdout_, STDOUT_FILENO);
            close(saved_stdout_);
            saved_stdout_ = -1;
        }
        if (fd_ >= 0) {
            close(fd_);
            fd_ = -1;
        }
        if (!path_.empty() && !finished_) {
            unlink(path_.c_str());
        }
    }

    int fd_ = -1;
    int saved_stdout_ = -1;
    bool finished_ = false;
    std::string path_;
    std::string output_;
};

} // namespace

std::vector<std::string> ScientificSimulationBackend::build_event_args(const GenulensConfig &config) const
{
    const auto &imf = config.model.imf;
    std::vector<std::string> args = {
        "genulens",
        "l", std::to_string(config.l),
        "b", std::to_string(config.b),
        "Nsimu", std::to_string(config.n_simu),
        "seed", std::to_string(config.seed),
        "tE", std::to_string(config.observed_tE), std::to_string(config.observed_tE_error),
        "M0", std::to_string(imf.m0),
        "M1", std::to_string(imf.m1),
        "M2", std::to_string(imf.m2),
        "M3", std::to_string(imf.m3),
        "Ml", std::to_string(imf.ml),
        "Mu", std::to_string(imf.mu),
        "alpha0", std::to_string(imf.alpha0),
        "alpha1", std::to_string(imf.alpha1),
        "alpha2", std::to_string(imf.alpha2),
        "alpha3", std::to_string(imf.alpha3),
        "alpha4", std::to_string(imf.alpha4),
        "VERBOSITY", "3",
    };

    for (std::size_t i = 1; i < config.raw_cli_args.size(); ++i) {
        const auto &arg = config.raw_cli_args[i];
        if (arg == "VERBOSITY" || arg == "l" || arg == "b" || arg == "Nsimu" || arg == "seed" ||
            arg == "tE" || arg == "M0" || arg == "M1" || arg == "M2" || arg == "M3" ||
            arg == "Ml" || arg == "Mu" || arg == "alpha0" || arg == "alpha1" ||
            arg == "alpha2" || arg == "alpha3" || arg == "alpha4") {
            ++i;
            if (arg == "tE") ++i;
            continue;
        }
        args.push_back(arg);
    }
    return args;
}

std::string ScientificSimulationBackend::run_cli_capture(const std::vector<std::string> &args) const
{
    std::vector<char *> argv;
    argv.reserve(args.size());
    for (const auto &arg : args) {
        argv.push_back(const_cast<char *>(arg.c_str()));
    }

    StdoutCapture capture;
    const int code = run_scientific_cli(static_cast<int>(argv.size()), argv.data());
    const auto output = capture.finish();
    if (code != 0) {
        throw std::runtime_error("genulens scientific backend returned non-zero status");
    }
    return output;
}

SimulationResult ScientificSimulationBackend::parse_verbosity3_events(const std::string &output, LikelihoodFunction likelihood) const
{
    SimulationResult result;
    std::istringstream lines(output);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream row(line);
        Event event;
        double pi_en = 0.0;
        double pi_ee = 0.0;
        double mu_sl = 0.0;
        double mu_sb = 0.0;
        double il = 0.0;
        double kl = 0.0;
        int remnant_flag = 0;

        if (!(row >> event.weight
                  >> event.lens_mass_msun
                  >> event.lens_distance_pc
                  >> event.source_distance_pc
                  >> event.tE
                  >> event.thetaE
                  >> event.piE
                  >> pi_en
                  >> pi_ee
                  >> event.mu_rel_masyr
                  >> mu_sl
                  >> mu_sb
                  >> il
                  >> kl
                  >> event.source_component
                  >> event.lens_component
                  >> remnant_flag)) {
            continue;
        }
        if (likelihood) {
            event.weight *= likelihood(event);
        }
        result.events.push_back(event);
    }
    return result;
}

SimulationResult ScientificSimulationBackend::simulate(const GenulensConfig &config, LikelihoodFunction likelihood) const
{
    return parse_verbosity3_events(run_cli_capture(build_event_args(config)), std::move(likelihood));
}

} // namespace genulens
