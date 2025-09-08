/**
 * @file
 * @brief Defines mathematical constants and functions for both CPU and CUDA execution.
 * @details This header provides a hardware-agnostic API for common mathematical operations
 * like logarithms, square roots, and random number generation. It uses preprocessor
 * directives to select the appropriate implementation (standard C++ or CUDA) at compile time.
 */
#ifndef MQI_MATH_HPP
#define MQI_MATH_HPP

#include <moqui/base/mqi_common.hpp>

#include <cmath>
#include <mutex>
#include <random>

namespace mqi
{

const float near_zero          = 1e-7f;      ///< A small value used to avoid floating-point comparison issues.
const float min_step           = 1e-3f;      ///< The minimum step size for particle transport.
const float geometry_tolerance = 1e-3f;      ///< A tolerance value for geometry intersection calculations.
const float m_inf              = -HUGE_VALF; ///< Negative infinity.
const float p_inf              = HUGE_VALF;  ///< Positive infinity.

/**
 * @brief Performs 1D linear interpolation.
 * @tparam T The floating-point type (e.g., float or double).
 * @param[in] x The point at which to interpolate.
 * @param[in] x0 The first x-coordinate.
 * @param[in] x1 The second x-coordinate.
 * @param[in] y0 The first y-coordinate, corresponding to x0.
 * @param[in] y1 The second y-coordinate, corresponding to x1.
 * @return The interpolated y-value at x.
 */
template<typename T>
CUDA_DEVICE inline T
intpl1d(T x, T x0, T x1, T y0, T y1) {
    return (x1 == x0) ? y0 : y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

/**
 * @brief Calculates the natural logarithm. Wrapper for `log` or `logf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The natural logarithm of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_ln(T s);

/**
 * @brief Calculates the square root. Wrapper for `sqrt` or `sqrtf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The square root of s.
 */
template<typename T>
CUDA_HOST_DEVICE T
mqi_sqrt(T s);

/**
 * @brief Calculates a number raised to a power. Wrapper for `pow` or `powf`.
 * @tparam T The floating-point type.
 * @param[in] s The base.
 * @param[in] p The exponent.
 * @return s raised to the power of p.
 */
template<typename T>
CUDA_DEVICE T
mqi_pow(T s, T p);

/**
 * @brief Calculates the base-e exponential. Wrapper for `exp` or `expf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return e raised to the power of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_exp(T s);

/**
 * @brief Calculates the arc cosine. Wrapper for `acos` or `acosf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The arc cosine of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_acos(T s);

/**
 * @brief Calculates the cosine. Wrapper for `cos` or `cosf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The cosine of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_cos(T s);

/**
 * @brief Calculates the sine. Wrapper for `sin` or `sinf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The sine of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_sin(T s);

/**
 * @brief Calculates the absolute value. Wrapper for `abs` or `fabs`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The absolute value of s.
 */
template<typename T>
CUDA_DEVICE T
mqi_abs(T s);

/**
 * @brief Rounds a number to the nearest integer. Wrapper for `round` or `roundf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The rounded value.
 */
template<typename T>
CUDA_HOST_DEVICE T
mqi_round(T s);

/**
 * @brief Calculates the floor of a number. Wrapper for `floor` or `floorf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The floor value.
 */
template<typename T>
CUDA_HOST_DEVICE T
mqi_floor(T s);

/**
 * @brief Calculates the ceiling of a number. Wrapper for `ceil` or `ceilf`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return The ceiling value.
 */
template<typename T>
CUDA_HOST_DEVICE T
mqi_ceil(T s);

/**
 * @brief Checks if a number is NaN (Not a Number). Wrapper for `isnan`.
 * @tparam T The floating-point type.
 * @param[in] s The value.
 * @return True if s is NaN, false otherwise.
 */
template<typename T>
CUDA_HOST_DEVICE bool
mqi_isnan(T s);

/**
 * @struct rnd_return
 * @brief A helper struct to define the return type for random number generators.
 * @tparam T The desired floating-point type (float or double).
 */
template<class T>
struct rnd_return {
    typedef T type;
};

/// @brief Specialization for float.
template<>
struct rnd_return<float> {
    typedef float type;
};

/// @brief Specialization for double.
template<>
struct rnd_return<double> {
    typedef double type;
};

/**
 * @brief Generates a normally distributed random number.
 * @tparam T The floating-point type of the distribution parameters.
 * @tparam S The type of the random number generator state/engine.
 * @param[in,out] rng A pointer to the random number generator.
 * @param[in] avg The mean of the distribution.
 * @param[in] sig The standard deviation of the distribution.
 * @return A random number from the specified normal distribution.
 */
template<class T, class S>
CUDA_DEVICE typename rnd_return<T>::type
mqi_normal(S* rng, T avg, T sig) {
    return T();
}

/**
 * @brief Generates a uniformly distributed random number in [0, 1).
 * @tparam T The floating-point type of the return value.
 * @tparam S The type of the random number generator state/engine.
 * @param[in,out] rng A pointer to the random number generator.
 * @return A random number from the uniform distribution.
 */
template<class T, class S>
CUDA_DEVICE typename rnd_return<T>::type
mqi_uniform(S* rng) {
    return T();
}

/**
 * @brief Generates an exponentially distributed random number.
 * @tparam T The floating-point type of the distribution parameters.
 * @tparam S The type of the random number generator state/engine.
 * @param[in,out] rng A pointer to the random number generator.
 * @param[in] avg The mean of the distribution (lambda = 1/avg).
 * @param[in] up The upper bound for the generated random number.
 * @return A random number from the specified exponential distribution.
 */
template<class T, class S>
CUDA_DEVICE typename rnd_return<T>::type
mqi_exponential(S* rng, T avg, T up) {
    return T();
}

#if defined(__CUDACC__)

///< specialization of template functions
//Nautual log
template<>
float
mqi_ln(float s) {
    return logf(s);
}
template<>
double
mqi_ln(double s) {
    return log(s);
}

///< sqrt
template<>
float
mqi_sqrt(float s) {
    return sqrtf(s);
}
template<>
double
mqi_sqrt(double s) {
    return sqrt(s);
}

///< power
template<>
float
mqi_pow(float s, float p) {
    return powf(s, p);
}
template<>
double
mqi_pow(double s, double p) {
    return pow(s, p);
}

///< exponential
template<>
float
mqi_exp(float s) {
    return expf(s);
}
template<>
double
mqi_exp(double s) {
    return exp(s);
}

///< acos
template<>
float
mqi_acos(float s) {
    return acosf(s);
}
template<>
double
mqi_acos(double s) {
    return acos(s);
}

///< cos
template<>
float
mqi_cos(float s) {
    return cosf(s);
}
template<>
double
mqi_cos(double s) {
    return cos(s);
}
///< sin
template<>
float
mqi_sin(float s) {
    return sinf(s);
}
template<>
double
mqi_sin(double s) {
    return sin(s);
}

template<>
float
mqi_abs(float s) {
    return abs(s);
}
template<>
double
mqi_abs(double s) {
    return abs(s);
}

template<>
float
mqi_round(float s) {
    return roundf(s);
}
template<>
double
mqi_round(double s) {
    return round(s);
}

template<>
float
mqi_floor(float s) {
    return floorf(s);
}
template<>
double
mqi_floor(double s) {
    return floor(s);
}

template<>
float
mqi_ceil(float s) {
    return ceilf(s);
}
template<>
double
mqi_ceil(double s) {
    return ceil(s);
}

template<>
bool
mqi_isnan(float s) {
    bool t = isnan(s);
    return t;
}
template<>
bool
mqi_isnan(double s) {
    bool t = isnan(s);
    return t;
}

//random number status per thread. each status is initialized by master seed.
//the parameters to be passed are currand_status for CUDA and random_engine for C++
//curand_status == mqi_rng ;
typedef curandState_t mqi_rng;

template<>
float
mqi_uniform<float>(mqi_rng* rng) {
    return curand_uniform(rng);
}

template<>   //template<class S = curandState>
double
mqi_uniform<double>(mqi_rng* rng) {
    return curand_uniform_double(rng);
}

template<>
float
mqi_normal<float>(mqi_rng* rng, float avg, float sig) {
    return curand_normal(rng) * sig + avg;
}

template<>
double
mqi_normal<double>(mqi_rng* rng, double avg, double sig) {
    return curand_normal_double(rng) * sig + avg;
}

template<>
float
mqi_exponential<float>(mqi_rng* rng, float avg, float up) {
    float x;
    do {
        x = -1.0 / avg * logf(1.0 - curand_uniform(rng));   //0, up
    } while (x > up || mqi::mqi_isnan(x));

    return x;
}

template<>
double
mqi_exponential<double>(mqi_rng* rng, double avg, double up) {
    double x;
    do {
        x = -1.0 / avg * log(1.0 - curand_uniform(rng));   //0, up
    } while (x > up || mqi::mqi_isnan(x));
    return x;
}

#else

//Natural log. C++ casts float to double. they have same implementation.
template<>
float
mqi_ln(float s) {
    return std::log(s);
}

template<>
double
mqi_ln(double s) {
    return std::log(s);
}

template<>
float
mqi_sqrt(float s) {
    return std::sqrt(s);
}

template<>
double
mqi_sqrt(double s) {
    return std::sqrt(s);
}

template<>
float
mqi_pow(float s, float p) {
    return std::pow(s, p);
}

template<>
double
mqi_pow(double s, double p) {
    return std::pow(s, p);
}

///< exponential
template<>
float
mqi_exp(float s) {
    return std::exp(s);
}
template<>
double
mqi_exp(double s) {
    return std::exp(s);
}

///< acos
template<>
float
mqi_acos(float s) {
    return std::acos(s);
}
template<>
double
mqi_acos(double s) {
    return std::acos(s);
}

template<>
float
mqi_cos(float s) {
    return std::cos(s);
}
template<>
double
mqi_cos(double s) {
    return std::cos(s);
}

template<>
float
mqi_abs(float s) {
    return std::abs(s);
}
template<>
double
mqi_abs(double s) {
    return std::abs(s);
}

///< round
template<>
float
mqi_round(float s) {
    return std::roundf(s);
}
template<>
double
mqi_round(double s) {
    return std::round(s);
}

////< floor

template<>
float
mqi_floor(float s) {
    return std::floor(s);
}
template<>
double
mqi_floor(double s) {
    return std::floor(s);
}

////< ceil
template<>
float
mqi_ceil(float s) {
    return std::ceil(s);
}
template<>
double
mqi_ceil(double s) {
    return std::ceil(s);
}

////< isnan
template<>
bool
mqi_isnan(float s) {
    return std::isnan(s);
}
template<>
bool
mqi_isnan(double s) {
    return std::isnan(s);
}

//compile: ERROR
typedef std::default_random_engine mqi_rng;
//typedef std::mt19937 mqi_rng;

template<>
float
mqi_uniform<float>(mqi_rng* rng) {
    std::uniform_real_distribution<float> dist;
    return dist(*rng);
}

template<>
double
mqi_uniform<double>(mqi_rng* rng) {
    std::uniform_real_distribution<double> dist;
    return dist(*rng);
}

template<>
float
mqi_normal<float>(mqi_rng* rng, float avg, float sig) {
    std::normal_distribution<float> dist(avg, sig);
    return dist(*rng);
}

template<>
double
mqi_normal<double>(mqi_rng* rng, double avg, double sig) {
    std::normal_distribution<double> dist(avg, sig);
    return dist(*rng);
}

template<>
float
mqi_exponential<float>(mqi_rng* rng, float avg, float up) {
    float                                x;
    std::exponential_distribution<float> dist(avg);
    x = dist(*rng);
    //    do {
    //        x = dist(*rng);
    //    } while (x > up || x <= 0);
    return x;
}

template<>
double
mqi_exponential<double>(mqi_rng* rng, double avg, double up) {
    double                                x;
    std::exponential_distribution<double> dist(avg);
    do {
        x = dist(*rng);
    } while (x > up || x <= 0);
    return x;
}

#endif

}   // namespace mqi

#endif
