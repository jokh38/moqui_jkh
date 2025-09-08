/// \file
/// \brief Defines error-checking functions and macros for CUDA and CPU operations.
///
/// This file provides utilities for checking and handling errors, particularly
/// those originating from CUDA API calls. It includes functions to check the last
/// CUDA error and a macro for wrapping CUDA calls to automatically check their return status.

#ifndef MQI_ERROR_CHECK_HPP
#define MQI_ERROR_CHECK_HPP

namespace mqi
{
/// \brief Checks the last error returned by the CUDA runtime.
///
/// This function retrieves the last error that occurred in the CUDA runtime.
/// If an error is found, it prints the error message, current GPU memory usage,
/// and exits the program. This is useful for debugging asynchronous CUDA kernels.
/// \param[in] msg A custom message to print alongside the error string.
inline void
check_cuda_last_error(const char* msg) {

#if defined(__CUDACC__)

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess) {
        size_t free, total;
        cudaMemGetInfo(&free, &total);
        printf(
          "current total %f MB free %f MB\n", total / (1024.0 * 1024.0), free / (1024.0 * 1024.0));
        printf("CUDA error: %s %s\n", msg, cudaGetErrorString(err));
        exit(-1);
    }

#endif
}

}   // namespace mqi

#if defined(__CUDACC__)
/// \def gpu_err_chk(ans)
/// \brief A macro to check the result of a CUDA API call.
///
/// This macro wraps a CUDA API call and passes its return code to the
/// `cuda_error_checker` function, along with the file and line number
/// where the call was made.
///
/// Example:
/// ```cpp
/// gpu_err_chk(cudaMalloc(&d_data, size));
/// ```
/// \param ans The CUDA API call to be checked.
#define gpu_err_chk(ans)                                                                           \
    { cuda_error_checker((ans), __FILE__, __LINE__); }

/// \brief Checks a CUDA error code and reports it if it's not `cudaSuccess`.
///
/// This function checks the given `cudaError_t` code. If it indicates an error,
/// it prints an error message to `stderr`, including the error string, file,
/// and line number.
/// \param[in] code The `cudaError_t` code to check.
/// \param[in] file The name of the file where the error occurred.
/// \param[in] line The line number where the error occurred.
/// \param[in] abort If true, the program will exit upon encountering an error.
inline void
cuda_error_checker(cudaError_t code, const char* file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

#endif

#endif