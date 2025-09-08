/**
 * @file
 * @brief Defines a struct containing fundamental physical constants and unit conversions.
 */
#ifndef MQI_PHYSICS_CONSTANTS_HPP
#define MQI_PHYSICS_CONSTANTS_HPP

#include <moqui/base/mqi_common.hpp>

namespace mqi
{

/**
 * @struct physics_constants
 * @brief A collection of fundamental physical constants and unit conversions.
 * @tparam R The floating-point type (e.g., float or double) used for the constants.
 * @details This struct provides a centralized and consistent source for physical constants
 * required in the simulation, such as particle masses, lengths, and energy units.
 */
template<typename R>
struct physics_constants {
    const R mm = 1.0;   ///< Default length unit (millimeter).
    const R cm = 10.0;  ///< Centimeter in terms of millimeters.
    const R cm3 =
        cm * cm * cm;   ///< Cubic centimeter in terms of cubic millimeters.
    const R mm3   = mm * mm * mm;   ///< Cubic millimeter.
    const R MeV   = 1.0;            ///< Default energy unit (Mega-electron Volt).
    const R eV    = 1e-6;           ///< Electron Volt in terms of MeV.
    const R Mp    = 938.272046 * MeV;   ///< Proton mass in MeV.
    const R Mp_sq = Mp * Mp;            ///< Square of the proton mass.
    const R Me    = 0.510998928 * MeV;  ///< Electron mass in MeV.
    const R Mo    = 14903.3460795634 * MeV;   ///< Oxygen atom mass in MeV.
    const R Mo_sq = Mo * Mo;                  ///< Square of the oxygen mass.
    const R MoMp  = Mo / Mp;                  ///< Ratio of Oxygen mass to Proton mass.
    const R MoMp_sq = MoMp * MoMp;            ///< Square of the Oxygen/Proton mass ratio.
    const R MeMp = Me / Mp;     ///< Ratio of Electron mass to Proton mass.
    const R MeMp_sq = MeMp * MeMp;   ///< Square of the Electron/Proton mass ratio.
    const R re =
        2.8179403262e-12 * mm;   ///< Classical electron radius in mm.
    const R re_sq            = re * re;   ///< Square of the classical electron radius.
    const R two_pi_re2_mc2   = 2.0 * M_PI * re_sq * Me;   ///< A pre-calculated constant: 2 * pi * r_e^2 * m_e * c^2.
    const R two_pi_re2_mc2_h2o =
        two_pi_re2_mc2 * 3.3428e+23 / cm3;   ///< The `two_pi_re2_mc2` constant scaled by the electron density of water.
    const R water_density = 1.0 / cm3;   ///< Density of water (1 g/cm^3) in g/mm^3.
    const R radiation_length_water =
        36.0863 * cm;   ///< Radiation length of water in mm.
    const R mev_to_joule = 1.60218e-13;   ///< Conversion factor from MeV to Joules.
};

}   // namespace mqi
#endif
