# Moqui C++

**A Monte Carlo simulation toolkit for radiotherapy**

## Description

Moqui C++ is a powerful and flexible toolkit for performing Monte Carlo simulations of radiation therapy treatments. It is designed to be highly performant, with support for GPU acceleration via CUDA, and can be used to simulate a variety of treatment modalities, including conventional radiotherapy (photons and electrons) and ion therapy (protons).

This project provides the low-level components for building complex simulations, including models for beamlines, apertures, and patient geometries. It is designed to be used as a library in other applications or through its own command-line interface for running simulations.

## Features

-   **High-performance Monte Carlo engine:** Optimized for speed with support for GPU acceleration using CUDA.
-   **Support for multiple modalities:** Can simulate photon, electron, and proton treatments.
-   **DICOM-RT compatibility:** Can interpret and use data from DICOM-RT plan and treatment record files.
-   **Flexible geometry and beam modeling:** Provides a set of classes for defining complex treatment machine geometries and beam sources.
-   **Extensible:** Designed to be modular, allowing for the addition of new physics models, geometries, and scoring options.

## Getting Started

### Prerequisites

-   **C++ Compiler**: A compiler that supports C++11 or later (e.g., GCC, Clang, MSVC).
-   **CMake**: Version 3.10 or later. You can download it from the [official CMake website](https://cmake.org/download/).
-   **NVIDIA CUDA Toolkit**: (Optional) Version 10.0 or later is required for GPU acceleration. You can download it from the [NVIDIA Developer website](https://developer.nvidia.com/cuda-toolkit).
-   **GDCM (Grassroots DICOM)**: A library for DICOM file handling. You can install it using a package manager.
    -   On Debian/Ubuntu: `sudo apt-get install libgdcm-tools`
    -   On macOS (using Homebrew): `brew install gdcm`
    -   Alternatively, download from the [GDCM SourceForge page](http://gdcm.sourceforge.net/).

### Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd moqui-cpp
    ```

2.  **Create a build directory:** It's best practice to build the project in a separate directory.
    ```bash
    mkdir build
    cd build
    ```

3.  **Run CMake to configure the project:**
    ```bash
    cmake ..
    ```
    If you have the CUDA Toolkit installed and want to enable GPU support, you can add the following option:
    ```bash
    cmake -DMOQUI_WITH_CUDA=ON ..
    ```

4.  **Build the project:**
    ```bash
    make
    ```
    After a successful build, the executables can be found in the `build/bin` directory.

## Usage

The project can be used as a library in your own C++ applications or through the provided command-line tools.

### As a Library

To use Moqui C++ as a library, you can include the necessary headers from the `base/` and other directories in your source files and link against the compiled library generated during the build process. The main entry point for simulations is typically through the `mqi::treatment_session` class.

### Command-Line Interface

The project includes several command-line tools for running simulations and tests. After building, these tools can be found in the `build/bin` directory.

A typical simulation is run using an input file that specifies the parameters for the simulation. For example:

```bash
./build/bin/executable_name -i /path/to/your/input.txt
```

An example `input.txt` file might look like this:

```
# Simulation Parameters
GPUID                   0
RandomSeed              12345
ParentDir               /path/to/patient/data
DicomDir                DICOM
logFilePath             Logs
OutputDir               /path/to/output/dir
OutputFormat            mhd
OverwriteResults        true
SimulationType          perBeam
BeamNumbers             1,2,3
ParticlesPerHistory     1000.0
SaveDoseToText          true
ScoreLet                true
```

#### Common Input Parameters

| Parameter           | Description                                                                                             |
| ------------------- | ------------------------------------------------------------------------------------------------------- |
| `GPUID`             | The ID of the GPU to use for the simulation.                                                            |
| `RandomSeed`        | The seed for the random number generator.                                                               |
| `ParentDir`         | The parent directory containing the patient data.                                                       |
| `DicomDir`          | The subdirectory within `ParentDir` that contains the DICOM files.                                      |
| `logFilePath`       | The path to the directory where log files will be saved.                                                |
| `OutputDir`         | The path to the directory where output files (e.g., dose grids) will be saved.                          |
| `OutputFormat`      | The format for the output files (e.g., `mhd`, `raw`).                                                     |
| `OverwriteResults`  | If `true`, existing output files will be overwritten.                                                   |
| `SimulationType`    | The type of simulation to run (`perBeam` or `perSpot`).                                                 |
| `BeamNumbers`       | A comma-separated list of beam numbers to simulate.                                                     |
| `ParticlesPerHistory` | The number of particles to simulate per history.                                                        |
| `SaveDoseToText`    | If `true`, saves the final dose grid as a text file in addition to the specified `OutputFormat`.        |
| `ScoreLet`          | If `true`, enables the calculation and scoring of Linear Energy Transfer (LET).                         |

For a full list of available options and parameters, please refer to the source code and the Doxygen documentation.

## Documentation

Comprehensive documentation for the Moqui C++ API can be generated using [Doxygen](https://www.doxygen.nl/). To generate the documentation, you will need to have Doxygen installed.

- On Debian/Ubuntu: `sudo apt-get install doxygen graphviz`
- On macOS (using Homebrew): `brew install doxygen`

Once installed, run the following command from the root of the repository:
```bash
doxygen Doxyfile
```

This will generate an HTML documentation in the `docs/html` directory. Open `docs/html/index.html` in your browser to view the documentation.

The `code_structure.md` file provides a high-level overview of the repository's structure, and `progress.md` tracks the documentation status of the files. The `todolist.txt` file outlines the remaining documentation work.

## Testing

A formal testing suite has not yet been established for this project. Contributions in this area, such as adding unit tests or integration tests, are highly welcome.

## Project Structure

The repository is organized into the following main directories:

-   `base/`: Contains the core classes and data structures for the Moqui toolkit. This includes fundamental components like particle tracks, vertices, 3D grids, and physics interaction models.
    -   `distributions/`: Classes for sampling from various particle phase-space distributions (e.g., uniform, Gaussian).
    -   `environments/`: Classes that define the simulation world, such as phantoms or patient-specific CT-based environments.
    -   `materials/`: Classes for material definitions and handling material properties.
    -   `scorers/`: Classes for scoring physical quantities like dose, LET, and track length.
-   `kernel_functions/`: Contains CUDA kernels and other GPU-specific code for high-performance computations, such as particle transport.
-   `treatment_machines/`: Contains definitions for specific radiotherapy treatment machine models, including their geometric components.

A detailed breakdown of the file structure can be found in `code_structure.md`.

## Contributing

Contributions to Moqui C++ are welcome! If you would like to contribute, please follow these steps:

1.  Fork the repository.
2.  Create a new branch for your feature or bug fix.
3.  Make your changes and commit them with a clear and descriptive message.
4.  Push your changes to your fork.
5.  Create a pull request to the main repository.

Please ensure that your code adheres to the existing coding style and that you have added appropriate documentation and tests.

### For AI Agents

If you are an AI agent working on this repository, please look for `AGENTS.md` files in the directory structure. These files may contain specific instructions or conventions to follow.

## License

This project is licensed under the [MIT License](LICENSE).
