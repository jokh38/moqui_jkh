#ifndef MQI_RANGESHIFTER_H
#define MQI_RANGESHIFTER_H

/// \file
///
/// RT-Ion geometry for rangeshifter

#include <moqui/base/mqi_geometry.hpp>

namespace mqi
{

/**
 * @class rangeshifter
 * @brief Represents a range shifter geometry used in radiotherapy.
 * @details This class defines the geometry of a range shifter, which can be either a rectangular box or a cylinder.
 * It inherits from the base `geometry` class and stores its specific dimensions and shape type.
 * It is assumed to consist of a single material.
 */
class rangeshifter : public geometry
{

public:
    const bool             is_rectangle;   ///< Flag indicating the shape of the volume. `true` for a rectangle, `false` for a cylinder.
    const mqi::vec3<float> volume;         ///< The dimensions of the volume. For a rectangle: (x, y, z). For a cylinder: (radius, theta, thickness).

    /**
     * @brief Constructs a rangeshifter object.
     * @param v The dimensions of the volume (x,y,z for rectangle; r, theta, thickness for cylinder).
     * @param p The position of the rangeshifter's center.
     * @param r The rotation matrix defining the orientation of the rangeshifter.
     * @param is_rect A boolean flag, `true` if the rangeshifter is rectangular (default), `false` if cylindrical.
     */
    rangeshifter(mqi::vec3<float>&   v,   ///< x,y,z or r, theta, thickness
                 mqi::vec3<float>&   p,   ///< position
                 mqi::mat3x3<float>& r,   ///< rotation matrix
                 bool                is_rect = true) :
        volume(v),
        is_rectangle(is_rect), geometry(p, r, mqi::geometry_type::RANGESHIFTER) {
        ;
    }

    /**
     * @brief Copy constructor.
     * @param rhs The rangeshifter object to copy.
     */
    rangeshifter(const mqi::rangeshifter& rhs) :
        volume(rhs.volume), is_rectangle(rhs.is_rectangle),
        geometry(rhs.pos, rhs.rot, rhs.geotype) {
        ;
    }

    /**
     * @brief Destructor.
     */
    ~rangeshifter() {
        ;
    }

    /**
     * @brief Assignment operator.
     * @param rhs The rangeshifter object to assign from.
     * @return A const reference to the assigned object.
     * @note This operator currently returns the rhs object directly and does not perform a proper assignment.
     */
    const rangeshifter&
    operator=(const mqi::rangeshifter& rhs) {
        return rhs;
    }
};

}   // namespace mqi

#endif