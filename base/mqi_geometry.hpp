/// \file
///
/// \brief Defines the base class for all geometric components in the simulation.
///
/// This file contains the definition of the `mqi::geometry` class, which serves
/// as the abstract base class for all geometric objects used in the Moqui
/// simulation environment, such as beamline components and patient anatomy.

#ifndef MQI_GEOMETRY_HPP
#define MQI_GEOMETRY_HPP

#include <array>
#include <map>
#include <string>
#include <vector>

#include <moqui/base/mqi_matrix.hpp>
#include <moqui/base/mqi_vec.hpp>

namespace mqi
{

/// \enum geometry_type
/// \brief Enumerates the different types of geometric components.
/// \note Many of these types are not yet fully supported.
typedef enum
{
    SNOUT,        ///< A component of the beam delivery system.
    RANGESHIFTER, ///< A range shifter.
    COMPENSATOR,  ///< A patient-specific compensator.
    BLOCK,        ///< A block aperture.
    BOLI,         ///< A bolus.
    WEDGE,        ///< A physical or dynamic wedge.
    TRANSFORM,    ///< A coordinate transformation.
    MLC,          ///< A multi-leaf collimator.
    PATIENT,      ///< The patient geometry.
    DOSEGRID,     ///< A dose calculation grid.
    UNKNOWN1,     ///< Unspecified geometry type 1.
    UNKNOWN2,     ///< Unspecified geometry type 2.
    UNKNOWN3,     ///< Unspecified geometry type 3.
    UNKNOWN4      ///< Unspecified geometry type 4.
} geometry_type;

/// \class geometry
/// \brief An abstract base class for all geometric objects in the simulation.
///
/// This class provides the fundamental properties of any geometric component,
/// including its position, orientation (rotation), and type. All specific
/// geometry classes should inherit from this class.
class geometry
{

public:
    ///< The position of the geometric object's origin.
    const mqi::vec3<float> pos;
    ///< The rotation matrix defining the object's orientation.
    const mqi::mat3x3<float> rot;
    ///< The type of the geometry, as defined by `geometry_type`.
    const geometry_type geotype;

    /// \brief Constructs a geometry object with a given position, rotation, and type.
    /// \param[in] p_xyz A reference to the position vector.
    /// \param[in] rot_xyz A reference to the 3x3 rotation matrix.
    /// \param[in] t The geometry type.
    geometry(mqi::vec3<float>& p_xyz, mqi::mat3x3<float>& rot_xyz, mqi::geometry_type t) :
        pos(p_xyz), rot(rot_xyz), geotype(t) {
        ;
    }

    /// \brief Constructs a geometry object with a given constant position, rotation, and type.
    /// \param[in] p_xyz A constant reference to the position vector.
    /// \param[in] rot_xyz A constant reference to the 3x3 rotation matrix.
    /// \param[in] t The geometry type.
    geometry(const mqi::vec3<float>&   p_xyz,
             const mqi::mat3x3<float>& rot_xyz,
             const mqi::geometry_type  t) :
        pos(p_xyz),
        rot(rot_xyz), geotype(t) {
        ;
    }

    /// \brief Copy constructor.
    /// \param[in] rhs The geometry object to copy.
    geometry(const geometry& rhs) : geotype(rhs.geotype), pos(rhs.pos), rot(rhs.rot) {
        ;
    }

    /// \brief Assignment operator.
    /// \param[in] rhs The geometry object to assign from.
    /// \return A constant reference to the assigned object.
    const geometry&
    operator=(const mqi::geometry& rhs) {
        return rhs;
    }

    /// \brief Virtual destructor.
    ~geometry() {
        ;
    }

    /// \brief A virtual method to print information about the geometry.
    ///
    /// Derived classes should override this method to provide specific details.
    virtual void
    dump() const {
        ;
    }
};

}   // namespace mqi
#endif
