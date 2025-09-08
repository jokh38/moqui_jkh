#ifndef MQI_COMMON_HPP
#define MQI_COMMON_HPP

/// \file mqi_common.hpp
///
/// \brief A header file containing common definitions, macros, and type aliases for the project.
///
/// This file includes CUDA-related headers and defines macros to handle code compilation
/// for both CPU and GPU, ensuring cross-platform compatibility. It also defines several
/// type aliases and enumerations used throughout the Moqui codebase.

#if defined(__CUDACC__)

#include <cublas.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cuda_fp16.h>
#include <curand.h>
#include <nvfunctional>
#include <stdio.h>

/// \def CUDA_HOST_DEVICE
/// \brief A macro that expands to `__host__ __device__` when compiling with NVCC,
/// allowing a function to be compiled for both CPU and GPU.
#define CUDA_HOST_DEVICE __host__ __device__
/// \def CUDA_HOST
/// \brief A macro that expands to `__host__` when compiling with NVCC.
#define CUDA_HOST __host__
/// \def CUDA_DEVICE
/// \brief A macro that expands to `__device__` when compiling with NVCC.
#define CUDA_DEVICE __device__
/// \def CUDA_GLOBAL
/// \brief A macro that expands to `__global__` when compiling with NVCC.
#define CUDA_GLOBAL __global__
/// \def CUDA_LOCAL
/// \brief A macro that expands to `__local__` when compiling with NVCC.
#define CUDA_LOCAL __local__
/// \def CUDA_SHARED
/// \brief A macro that expands to `__shared__` when compiling with NVCC.
#define CUDA_SHARED __shared__
/// \def CUDA_CONSTANT
/// \brief A macro that expands to `__constant__` when compiling with NVCC.
#define CUDA_CONSTANT __constant__

#else

#define CUDA_HOST_DEVICE
#define CUDA_HOST
#define CUDA_DEVICE
#define CUDA_GLOBAL
#define CUDA_LOCAL
#define CUDA_SHARED
#define CUDA_CONSTANT
#endif

#include <cmath>
#include <cstdint>
#include <limits>

namespace mqi
{
/// \typedef phsp_t
/// \brief Type alias for phase-space variables, typically float.
typedef float phsp_t;

/// \typedef cnb_t
/// \brief Type alias for cumulative number of something (e.g., histories), typically uint64_t.
typedef uint64_t cnb_t;
/// \typedef ijk_t
/// \brief Type alias for grid indices (i, j, k), typically int32_t.
typedef int32_t ijk_t;

#if defined(__CUDACC__)
/// \typedef density_t
/// \brief Type alias for density values, typically float. On CUDA devices, this could be `__half`.
typedef float density_t;
#else
typedef float density_t;
#endif

/// \brief The maximum number of blocks in a CUDA grid (limited to 2^16-1).
const uint16_t block_limit = 65535;
/// \brief The maximum number of threads in a CUDA block.
const uint16_t thread_limit = 512;

/// \typedef key_t
/// \brief Type alias for keys used in hash tables or maps, typically uint32_t.
typedef uint32_t key_t;
/// \brief A constant representing an empty or invalid key.
const key_t empty_pair = 0xffffffff;

/// \enum cell_side
/// \brief Enumerates the six possible faces of a voxel that can be intersected.
typedef enum
{
    XM = 0, ///< The face on the minimum-x side.
    XP = 1, ///< The face on the maximum-x side.
    YM = 2, ///< The face on the minimum-y side.
    YP = 3, ///< The face on the maximum-y side.
    ZM = 4, ///< The face on the minimum-z side.
    ZP = 5, ///< The face on the maximum-z side.
    NONE_XYZ_PLANE = 6 ///< Represents no intersection or an invalid side.
} cell_side;

/// \enum aperture_type_t
/// \brief Enumerates the types of apertures.
typedef enum
{
    MASK   = 0, ///< A mask-type aperture.
    VOLUME = 1  ///< A volume-type aperture.
} aperture_type_t;

/// \enum sim_type_t
/// \brief Enumerates the different simulation modes.
typedef enum
{
    PER_BEAM    = 0, ///< Simulation is run on a per-beam basis.
    PER_SPOT    = 1, ///< Simulation is run on a per-spot basis.
    PER_PATIENT = 2  ///< Simulation is run for the entire patient.
} sim_type_t;

/// \enum transport_type
/// \brief Enumerates different particle transport conditions.
typedef enum
{
    APERTURE_CLOSE = 1, ///< Particle is outside the aperture opening.
    APERTURE_OPEN  = 2, ///< Particle is inside the aperture opening.
    NORMAL_PHYSICS = 3  ///< Standard physics transport, not interacting with an aperture.
} transport_type;

/// \brief A global constant for the maximum step size in the simulation.
const float max_step_global = 1.0;

}   // namespace mqi

#endif
