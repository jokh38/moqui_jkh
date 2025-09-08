#ifndef MQI_BEAMLET_HPP
#define MQI_BEAMLET_HPP

/// \file mqi_beamlet.hpp
///
/// \brief Defines a beamlet, a fundamental component of a beam model.
///
/// A beamlet represents a small, idealized component of a radiation beam, such as a
/// pencil beam in intensity-modulated proton therapy (IMPT). It is defined by a
/// collection of probability distributions that model the phase-space of particles
/// (energy, position, and direction).

#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <tuple>

#include <moqui/base/mqi_coordinate_transform.hpp>
#include <moqui/base/mqi_distributions.hpp>
#include <moqui/base/mqi_vertex.hpp>

namespace mqi
{

/// \class beamlet
/// \brief Represents a single beamlet within a treatment field.
///
/// A beamlet can model a small part of a larger field (like a pencil beam in IMPT)
/// or an entire field (like a double-scattered proton field). It combines an energy
/// distribution and a 6D fluence distribution (x, x', y, y', z, z') to generate
/// initial particle states (vertices). It also includes a coordinate transformation
/// to map the locally generated particles into the patient or treatment coordinate system.
///
/// \tparam T The data type for numerical values (e.g., float, double).
/// \note This class is designed to be usable on both CPU and GPU (via CUDA).
template<typename T>
class beamlet
{
protected:
    /// \brief A pointer to the 1D probability distribution function for particle energy.
    mqi::pdf_Md<T, 1>* energy = nullptr;

    /// \brief A pointer to the 6D probability distribution function for particle fluence.
    /// This distribution samples the phase-space variables (x, y, z, x', y', z').
    mqi::pdf_Md<T, 6>* fluence = nullptr;

    /// \brief The coordinate transformation to map from the local beamlet frame
    /// to the global patient or treatment coordinate system.
    coordinate_transform<T> p_coord;

public:
    /// \brief Constructs a beamlet from energy and fluence distributions.
    ///
    /// \param e A pointer to a 1D PDF for energy.
    /// \param f A pointer to a 6D PDF for fluence.
    CUDA_HOST_DEVICE
    beamlet(mqi::pdf_Md<T, 1>* e, mqi::pdf_Md<T, 6>* f) : energy(e), fluence(f) {
        ;
    }

    /// \brief Default constructor.
    CUDA_HOST_DEVICE
    beamlet() {
        ;
    }

    /// \brief Copy constructor.
    ///
    /// Creates a new beamlet as a copy of an existing one.
    /// \param rhs The beamlet object to copy.
    CUDA_HOST_DEVICE
    beamlet(const beamlet<T>& rhs) {
        energy  = rhs.energy;
        fluence = rhs.fluence;
        p_coord = rhs.p_coord;
    }

    /// \brief Sets the coordinate transformation for the beamlet.
    ///
    /// This function allows setting the rotation and translation that will be applied
    /// to the generated particles to place them in the correct global coordinate system.
    ///
    /// \param p A `coordinate_transform` object containing the desired translation and rotation.
    CUDA_HOST_DEVICE
    void
    set_coordinate_transform(coordinate_transform<T> p) {
        p_coord = p;
    }

    /// \brief Samples a particle vertex from the beamlet's distributions.
    ///
    /// This operator uses the energy and fluence distributions to generate a single
    /// particle's kinematic properties (energy, position, direction), applies the
    /// coordinate transformation, and returns the result as a vertex.
    ///
    /// \param rng A pointer to a random number engine to be used for sampling.
    /// \return A `mqi::vertex_t<T>` object representing the generated particle.
    virtual mqi::vertex_t<T>
    operator()(std::default_random_engine* rng) {
        std::array<T, 6> phsp = (*fluence)(rng);
        mqi::vec3<T>     pos(phsp[0], phsp[1], phsp[2]);
        mqi::vec3<T>     dir(phsp[3], phsp[4], phsp[5]);
        mqi::vertex_t<T> vtx;
        vtx.ke  = (*energy)(rng)[0];
        vtx.pos = p_coord.rotation * pos + p_coord.translation;
        vtx.dir = p_coord.rotation * dir;
        return vtx;
    };
};

}   // namespace mqi

#endif
