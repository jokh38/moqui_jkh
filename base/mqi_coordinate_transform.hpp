/// \file
///
/// \brief Defines a coordinate transformation class for mapping points between coordinate systems.
///
/// This file contains the definition of the `coordinate_transform` class, which is
/// used to handle transformations (rotation and translation) between different
/// coordinate systems, such as from a beam source to IEC or DICOM coordinates.

#ifndef MQI_COORDINATE_TRANSFORM_HPP
#define MQI_COORDINATE_TRANSFORM_HPP

#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <tuple>

#include <moqui/base/mqi_matrix.hpp>
#include <moqui/base/mqi_vec.hpp>

namespace mqi
{

/// \class coordinate_transform
/// \brief Manages coordinate transformations using rotation and translation.
/// \tparam T The floating-point type for calculations.
///
/// This class handles the transformation of coordinates between different systems,
/// applying a series of rotations (collimator, gantry, couch) and a translation.
/// The default units are millimeters for length and degrees for angles. Rotations
/// are applied counter-clockwise (CCW). The order of rotation is:
/// collimator -> gantry -> couch -> iec2dicom.
template<typename T>
class coordinate_transform
{
public:
    ///< A constant to convert degrees to radians.
    const float deg2rad = M_PI / 180.0;

    ///< The translation vector for the transformation.
    vec3<T> translation;

    ///< An array of rotation angles [collimator, gantry, couch, iec].
    std::array<T, 4> angles = { 0.0, 0.0, 0.0, 0.0 };

    ///< Rotation matrix for the collimator.
    mat3x3<T> collimator;
    ///< Rotation matrix for the gantry.
    mat3x3<T> gantry;
    ///< Rotation matrix for the patient support (couch).
    mat3x3<T> patient_support;
    ///< Rotation matrix for IEC to DICOM transformation.
    mat3x3<T> iec2dicom;

    ///< The final combined rotation matrix.
    mat3x3<T> rotation;

    /// \brief Constructs a coordinate_transform object with specified angles and position.
    /// \param[in] ang A reference to an array of 4 angles (collimator, gantry, couch, iec).
    /// \param[in] pos A reference to the translation vector (isocenter).
    CUDA_HOST_DEVICE
    coordinate_transform(std::array<T, 4>& ang, vec3<T>& pos) : angles(ang), translation(pos) {
        collimator      = mat3x3<T>(0, 0, angles[0] * deg2rad);
        gantry          = mat3x3<T>(0, angles[1] * deg2rad, 0);
        patient_support = mat3x3<T>(0, 0, angles[2] * deg2rad);
        iec2dicom       = mat3x3<T>(angles[3] * deg2rad, 0, 0);
        rotation        = iec2dicom * patient_support * gantry * collimator;
    }

    /// \brief Constructs a coordinate_transform object with constant angles and position.
    /// \param[in] ang A constant reference to an array of 4 angles.
    /// \param[in] pos A constant reference to the translation vector.
    CUDA_HOST_DEVICE
    coordinate_transform(const std::array<T, 4>& ang, const vec3<T>& pos) :
        angles(ang), translation(pos) {
        collimator      = mat3x3<T>(0, 0, angles[0] * deg2rad);
        gantry          = mat3x3<T>(0, angles[1] * deg2rad, 0);
        patient_support = mat3x3<T>(0, 0, angles[2] * deg2rad);
        iec2dicom       = mat3x3<T>(angles[3] * deg2rad, 0, 0);
        rotation        = iec2dicom * patient_support * gantry * collimator;
    }

    /// \brief Copy constructor.
    /// \param[in] ref A constant reference to the coordinate_transform object to copy.
    CUDA_HOST_DEVICE
    coordinate_transform(const coordinate_transform<T>& ref) {
        angles          = ref.angles;
        collimator      = ref.collimator;
        gantry          = ref.gantry;
        patient_support = ref.patient_support;
        iec2dicom       = ref.iec2dicom;
        rotation        = iec2dicom * patient_support * gantry * collimator;
        translation     = ref.translation;
    }

    /// \brief Default constructor.
    CUDA_HOST_DEVICE
    coordinate_transform() {
        ;
    }

    /// \brief Assignment operator.
    /// \param[in] ref A constant reference to the coordinate_transform object to assign.
    /// \return A reference to the updated coordinate_transform object.
    CUDA_HOST_DEVICE
    coordinate_transform<T>&
    operator=(const coordinate_transform<T>& ref) {
        collimator      = ref.collimator;
        gantry          = ref.gantry;
        patient_support = ref.patient_support;
        iec2dicom       = ref.iec2dicom;
        rotation        = iec2dicom * patient_support * gantry * collimator;
        translation     = ref.translation;
        angles          = ref.angles;
        return *this;
    }

    /// \brief Prints the transformation matrices and vectors to the console.
    ///
    /// This method is available on the host and is useful for debugging.
    CUDA_HOST
    void
    dump() {
        std::cout << "--- coordinate transform---" << std::endl;
        std::cout << "    translation ---" << std::endl;
        translation.dump();
        std::cout << "    rotation ---" << std::endl;
        rotation.dump();
        std::cout << "    gantry ---" << std::endl;
        gantry.dump();
        std::cout << "    patient_support ---" << std::endl;
        patient_support.dump();
        std.cout << "    collimator ---" << std::endl;
        collimator.dump();
        std::cout << "    IEC2DICOM ---" << std::endl;
        iec2dicom.dump();
    }
};

}   // namespace mqi

#endif
