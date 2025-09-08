#ifndef MQI_CLI_HPP
#define MQI_CLI_HPP

/// \file mqi_cli.hpp
///
/// \brief Defines a command-line interface helper class.
///
/// This file provides the `cli` class, which is designed to parse and manage
/// command-line arguments for unit tests and other applications within the project.

#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

#include <sys/stat.h>
#include <cstdlib>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

namespace mqi
{

/// \class cli
/// \brief A simple command-line interface for unit tests and applications.
///
/// This class provides a basic framework for parsing command-line arguments.
/// It predefines a set of expected options and stores their corresponding values.
/// Its main purposes are:
/// 1. Reading RT-Ion plan files.
/// 2. Creating geometries, patients, dose grids, and beamline components.
/// 3. Creating beam sources for a machine model.
class cli
{

protected:
    /// \brief A map to store the command-line options and their arguments.
    /// The key is the option name (e.g., "--dicom_path"), and the value is a
    /// vector of string arguments provided for that option.
    std::map<const std::string, std::vector<std::string>> parameters;

public:
    /// \brief Constructs the `cli` object and initializes the list of expected parameters.
    cli() {
        parameters = std::map<const std::string, std::vector<std::string>>({
          { "--dicom_path", {} },
          { "--bname", {} },
          { "--bnumber", {} },
          { "--spots", {} },
          { "--beamlets", {} },
          { "--pph", {} },
          { "--sid", {} },
          { "--output_prefix", {} },
          { "--nhistory", {} },
          { "--pxyz", {} },
          { "--source_energy", {} },
          { "--energy_variance", {} },
          { "--rxyz", {} },
          { "--lxyz", {} },
          { "--nxyz", {} },
          { "--spot_position", {} },
          { "--spot_size", {} },
          { "--spot_angles", {} },
          { "--spot_energy", {} },
          { "--histories", {} },
          { "--threads", {} },
          { "--score_variance", {} },
          { "--gpu_id", {} },
          { "--output_format", {} },
          { "--random_seed", {} },
          { "--phantom_path", {} }
        });
    }

    /// \brief Destroys the `cli` object.
    ~cli() {
        ;
    }

    /// \brief Parses the command-line arguments and populates the parameters map.
    ///
    /// This function iterates through the command-line arguments (`argv`). When it
    /// encounters a known option (one that exists as a key in the `parameters` map),
    /// it consumes all subsequent arguments until it finds the next option (a string
    /// starting with "--"). These consumed arguments are stored as a vector of strings
    /// associated with the known option key.
    ///
    /// \param argc The number of command-line arguments.
    /// \param argv An array of C-style strings containing the command-line arguments.
    void
    read(int argc, char** argv) {
        std::cout << "# of arguments: " << argc << std::endl;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            auto it = parameters.find(arg);

            if (it != parameters.end()) {
                // This is a known option, so consume its arguments
                int j = i + 1;
                while (j < argc && std::string(argv[j]).rfind("--", 0) != 0) {
                    it->second.push_back(argv[j]);
                    j++;
                }
                i = j - 1; // The outer loop will increment i to j
            } else {
                std::cerr << "Warning: Unknown option '" << arg << "' ignored." << std::endl;
            }
        }

        // Print out the parsed parameters for verification
        for (const auto& pair : parameters) {
            if (!pair.second.empty()) {
                std::cout << pair.first << " : ";
                for (const auto& param : pair.second) {
                    std::cout << param << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    /// \brief Prints a help message describing the available command-line options.
    ///
    /// \param s The name of the executable (`argv[0]`).
    virtual void
    print_help(char* s) {
        std::cout << "Usage:   " << s << " [-option] [argument]" << std::endl;
        std::cout << "options:  " << std::endl;
        for (const auto& pair : parameters) {
            std::cout << "          " << pair.first << std::endl;
        }
    }

    /// \brief Provides access to the arguments for a specific option.
    ///
    // \param t The name of the command-line option (e.g., "--dicom_path").
    /// \return A vector of strings containing the arguments for the specified option.
    ///         Returns an empty vector if the option was not provided.
    const std::vector<std::string>
    operator[](const std::string& t) {
        return parameters[t];
    }
};
}   // namespace mqi

#endif
