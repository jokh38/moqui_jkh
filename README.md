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

-   A C++ compiler that supports C++11 or later.
-   CMake (version 3.10 or later).
-   (Optional) NVIDIA CUDA Toolkit (version 10.0 or later) for GPU acceleration.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd moqui-cpp
    ```

2.  **Create a build directory:**
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

## Usage

The project can be used as a library in your own C++ applications or through the provided command-line tools.

### As a Library

To use Moqui C++ as a library, you can include the necessary headers in your source files and link against the compiled library.

### Command-Line Interface

The project includes several command-line tools for running simulations and tests. These tools can be found in the `build/bin` directory after building the project.

For example, to run a simulation using a DICOM-RT Ion Plan, you might use a command like this:
```bash
./path/to/executable --dicom_path /path/to/dicom --bname "MyBeam" --nhistory 1000000
```
For a full list of available options, you can run the executable with the `--help` flag.

## Project Structure

The repository is organized into the following main directories:

-   `base/`: Contains the core classes and data structures for the Moqui toolkit, such as geometries, beam models, and physics interactions.
-   `kernel_functions/`: Contains CUDA kernels and other GPU-related code.
-   `treatment_machines/`: Contains definitions for specific treatment machine models.

A detailed breakdown of the file structure can be found in `code_structure.md`.

## Contributing

Contributions to Moqui C++ are welcome! If you would like to contribute, please follow these steps:

1.  Fork the repository.
2.  Create a new branch for your feature or bug fix.
3.  Make your changes and commit them with a clear and descriptive message.
4.  Push your changes to your fork.
5.  Create a pull request to the main repository.

Please ensure that your code adheres to the existing coding style and that you have added appropriate documentation and tests.

## License

This project is licensed under the [MIT License](LICENSE).
