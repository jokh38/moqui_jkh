#ifndef MQI_REL_QUANTITIES_HPP
#define MQI_REL_QUANTITIES_HPP

#include <moqui/base/mqi_common.hpp>
#include <moqui/base/mqi_math.hpp>

namespace mqi
{

/**
 * @class relativistic_quantities
 * @brief A struct to calculate and store essential relativistic quantities for a particle.
 * @details This class takes a particle's kinetic energy and rest mass, and from these, it computes
 * several derived quantities like the Lorentz factor (gamma), beta, total energy, and maximum
 * energy transfer to an electron. This pre-calculation is an optimization to avoid repeated computations
 * in physics models.
 * @tparam R The floating-point type for calculations (e.g., float, double).
 */
template<typename R>
class relativistic_quantities
{
public:
    R beta_sq;    ///< The square of the particle's velocity relative to the speed of light (v^2/c^2).
    R gamma_sq;   ///< The square of the Lorentz factor.
    R gamma;      ///< The Lorentz factor (1 / sqrt(1 - beta^2)).
    R Et;         ///< The total energy of the particle (kinetic energy + rest mass energy).
    R Et_sq;      ///< The square of the total energy.
    R Ek;         ///< The kinetic energy of the particle in MeV.
    R mc2;        ///< The rest mass energy of the particle in MeV (m0*c^2).
    R tau;        ///< The ratio of kinetic energy to rest mass energy (Ek / mc^2).
    R Te_max;     ///< The maximum kinetic energy that can be transferred to a stationary electron in a single collision.

    /**
     * @brief Constructs and computes the relativistic quantities.
     * @details Initializes all member variables based on the provided kinetic energy and rest mass.
     * Note that the calculation currently assumes the incident particle is a proton for calculating `Et`, `gamma`, and `tau`.
     * The rest mass of the electron is also hardcoded for the `Te_max` calculation.
     * @param kinetic_energy The kinetic energy of the particle in MeV.
     * @param rest_mass_MeV The rest mass of the particle in MeV.
     */
    CUDA_HOST_DEVICE
    relativistic_quantities(R kinetic_energy, R rest_mass_MeV) :
        Ek(kinetic_energy), mc2(rest_mass_MeV) {
        const R Mp      = 938.272046;   // Mp/eV = proton mass in eV
        const R Me      = 0.510998928;
        Et              = Ek + Mp;
        Et_sq           = Et * Et;
        const R MeMp    = Me / Mp;       // Me/Mp
        const R MeMp_sq = MeMp * MeMp;   // (Me/Mp)^2
        gamma           = Et / Mp;
        gamma_sq        = gamma * gamma;
        beta_sq         = 1.0 - 1.0 / gamma_sq;

        Te_max = (2.0 * Me * beta_sq * gamma_sq);
        Te_max /= (1.0 + 2.0 * gamma * MeMp + MeMp_sq);
        tau = Ek / Mp;
    }

    /**
     * @brief Destructor.
     */
    CUDA_HOST_DEVICE
    ~relativistic_quantities() {
        ;
    }

    /**
     * @brief Calculates the relativistic momentum of the particle.
     * @return The momentum of the particle in units of MeV/c.
     */
    CUDA_HOST_DEVICE
    R
    momentum() {
        return mqi::mqi_sqrt(Et * Et - mc2 * mc2);
    }
};

}   // namespace mqi
#endif
