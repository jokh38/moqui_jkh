#ifndef MQI_INTERACTION_HPP
#define MQI_INTERACTION_HPP
#include <random>

#include <moqui/base/mqi_material.hpp>
#include <moqui/base/mqi_math.hpp>
#include <moqui/base/mqi_physics_constants.hpp>
#include <moqui/base/mqi_relativistic_quantities.hpp>
#include <moqui/base/mqi_track.hpp>
#include <moqui/base/mqi_track_stack.hpp>

namespace mqi
{

/**
 * @class interaction
 * @brief A pure virtual class representing the interaction between a particle and a material.
 * @tparam R The floating-point type (e.g., float or double).
 * @tparam P The particle type.
 * @details This class serves as an interface for different physics interaction models. It defines the common methods that any interaction process should implement, such as calculating cross-sections, sampling step lengths, and updating particle tracks.
 */
template<typename R, mqi::particle_t P>
class interaction
{
public:
    const physics_constants<R> units;   ///< Physics constants with appropriate units.
#ifdef __PHYSICS_DEBUG__
    R T_cut = 0.08511 * units.MeV;   ///< Kinetic energy cut for secondary particles in debug mode.
#else
    R T_cut = 0.0815 * units.MeV;   ///< Kinetic energy cut for secondary particles.
#endif
    const R max_step = 0.01 * units.cm;   ///< Maximum step size used to generate dE/dx and cross-section data.
    R       Tp_cut   = 0.5 * units.MeV;   ///< Kinetic energy cut for protons.
    const mqi::vec3<R>
        dir_z;   ///< Momentum direction of the incident particle on the scattering plane.

public:
    /**
     * @brief Default constructor.
     */
    CUDA_HOST_DEVICE
    interaction() : dir_z(0, 0, -1) {
        ;
    }

    /**
     * @brief Destructor.
     */
    CUDA_HOST_DEVICE
    ~interaction() {
        ;
    }

    /**
     * @brief Samples a step length for the particle.
     * @param[in] rel Relativistic quantities of the particle.
     * @param[in] mat The material through which the particle is traveling.
     * @param[in,out] rng A pointer to the random number generator.
     * @return The sampled step length in cm.
     * @details The step length is sampled from an exponential distribution based on the material's properties and the interaction cross-section.
     */
    CUDA_HOST_DEVICE
    virtual R
    sample_step_length(const relativistic_quantities<R>& rel,
                       const material_t<R>&              mat,
                       mqi_rng*                          rng) {
        R cs   = mat.rho_mass * this->cross_section(rel, mat);
        R mfp  = (cs == 0.0) ? max_step : 1.0 / cs;
        R prob = mqi_uniform<R>(rng);
        return -1.0 * mfp * mqi_ln(prob);
    }

    /**
     * @brief Samples a step length for the particle given a cross-section.
     * @param[in] cs The total cross-section.
     * @param[in,out] rng A pointer to the random number generator.
     * @return The sampled step length in cm.
     * @details The step length is sampled from an exponential distribution. This is an overloaded version of sample_step_length.
     */
    CUDA_HOST_DEVICE
    virtual R
    sample_step_length(const R cs, mqi_rng* rng) {
        R mfp  = (cs == 0.0) ? max_step : 1.0 / cs;
        R prob = mqi_uniform<R>(rng);
        return -1.0 * mfp * mqi_ln(prob);
    }

    /**
     * @brief Calculates the cross-section for the interaction.
     * @param[in] rel Relativistic quantities of the particle.
     * @param[in] mat The material.
     * @return The cross-section in cm^2.
     * @details This is a pure virtual function that must be implemented by derived classes.
     */
    CUDA_HOST_DEVICE
    virtual R
    cross_section(const relativistic_quantities<R>& rel, const material_t<R>& mat) = 0;

    /**
     * @brief Updates the particle track during a step (continuous effects).
     * @param[in,out] trk The particle track to be updated.
     * @param[in,out] stk The secondary particle stack.
     * @param[in,out] rng A pointer to the random number generator.
     * @param[in] len The step length.
     * @param[in,out] mat The material.
     * @details This pure virtual function handles continuous processes like energy loss along a step.
     */
    CUDA_HOST_DEVICE
    virtual void
    along_step(track_t<R>&       trk,
               track_stack_t<R>& stk,
               mqi_rng*          rng,
               const R           len,
               material_t<R>&    mat) = 0;

    /**
     * @brief Updates the particle track at the end of a step (discrete effects).
     * @param[in,out] trk The particle track to be updated.
     * @param[in,out] stk The secondary particle stack.
     * @param[in,out] rng A pointer to the random number generator.
     * @param[in] len The step length.
     * @param[in,out] mat The material.
     * @param[in] score_local_deposit A flag to indicate whether to score energy locally.
     * @details This pure virtual function handles discrete events at the end of a step, such as creating secondary particles.
     */
    CUDA_HOST_DEVICE
    virtual void
    post_step(track_t<R>&       trk,
              track_stack_t<R>& stk,
              mqi_rng*          rng,
              const R           len,
              material_t<R>&    mat,
              bool              score_local_deposit) = 0;
};

}   // namespace mqi

#endif
