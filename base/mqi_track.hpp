#ifndef MQI_TRACK_HPP
#define MQI_TRACK_HPP

#include <moqui/base/mqi_node.hpp>
#include <moqui/base/mqi_vertex.hpp>

namespace mqi
{

/**
 * @enum particle_t
 * @brief Defines the types of particles that can be simulated.
 */
typedef enum
{
    PHOTON   = 0, ///< A photon particle.
    ELECTRON = 1, ///< An electron particle.
    PROTON   = 2, ///< A proton particle.
    NEUTRON  = 3  ///< A neutron particle.
} particle_t;

/**
 * @enum process_t
 * @brief Defines the physics or geometry process that limits a particle's transport step.
 */
typedef enum
{
    BEGIN    = 0, ///< Initial state at the beginning of a track.
    MAX_STEP = 1, ///< Step was limited by the maximum allowed step size.
    BOUNDARY = 2, ///< Step was limited by a geometry boundary crossing.
    CSDA     = 3, ///< Step was limited by the continuous slowing down approximation range.
    D_ION    = 4, ///< Step was limited by a delta-electron ionization event.
    PP_E     = 5, ///< Step was limited by a proton-proton elastic scattering event.
    PO_E     = 6, ///< Step was limited by a proton-oxygen elastic scattering event.
    PO_I     = 7, ///< Step was limited by a proton-oxygen inelastic scattering event.
    KILLED0  = 8, ///< Particle was killed for a generic reason.
    KILLED1  = 9  ///< Particle was killed for another generic reason.
} process_t;

/**
 * @enum status_t
 * @brief Defines the current lifecycle status of a particle track.
 */
typedef enum
{
    CREATED = 0, ///< The track has been created but not yet transported.
    STOPPED = 3  ///< The track has been stopped due to a physics process or energy cut-off and will no longer be transported.
} status_t;

/**
 * @class track_t
 * @brief Represents a particle's state and its path through the simulation.
 * @details This class is a fundamental data structure that encapsulates all information about a particle
 * during its simulation, including its type, status, energy, position, direction, and the physics process
 * that governed its most recent step. It represents a single step of the particle's journey, from a
 * pre-step vertex (`vtx0`) to a post-step vertex (`vtx1`).
 * @tparam R The floating-point type for coordinate and energy values (e.g., float, double).
 */
template<typename R>
class track_t
{
public:
    uint32_t  scorer_column; ///< An ID, often for a beamlet or source, used to create influence matrices (e.g., Dij).
    status_t  status;        ///< The current lifecycle status of the particle (e.g., CREATED, STOPPED).
    process_t process;       ///< The physics or geometry process that limited the last step.
    bool      primary;       ///< Flag indicating if this is a primary particle (`true`) or a secondary (`false`).
    particle_t particle;     ///< The type of the particle (e.g., PROTON, ELECTRON).

    vertex_t<R> vtx0; ///< The pre-step vertex (start of the current step).
    vertex_t<R> vtx1; ///< The post-step vertex (end of the current step).

    R dE       = 0.0;     ///< The total energy deposited along the current step.
    R local_dE = 0.0;     ///< The energy deposited locally (e.g., from delta electrons below the tracking cut-off).
    node_t<R>* c_node = nullptr; ///< A pointer to the current geometry node (volume) the track is in.

    intersect_t<R> its; ///< Stores information about the next geometry intersection.

    vec3<R> ref_vector = vec3<R>(0, 0, 1); ///< A reference vector (0,0,1) used as a basis for rotations.

    /**
     * @brief Default constructor. Initializes a primary particle at the beginning of its life.
     */
    CUDA_HOST_DEVICE
    track_t() :
        status(CREATED), process(BEGIN), primary(true), dE(0), scorer_column(0), local_dE(0) {
        ;
    }

    /**
     * @brief Constructs a track starting from a given vertex.
     * @param v The initial vertex for the track. Both `vtx0` and `vtx1` are set to this value.
     */
    CUDA_HOST_DEVICE
    track_t(const vertex_t<R>& v) :
        status(CREATED), process(BEGIN), primary(true), dE(0), scorer_column(0), local_dE(0) {
        vtx0 = v;
        vtx1 = v;
    }

    /**
     * @brief Constructs a track with fully specified initial state.
     */
    CUDA_HOST_DEVICE
    track_t(status_t    s,
            process_t   p,
            bool        is_p,
            particle_t  t,
            vertex_t<R> v0,
            vertex_t<R> v1,
            const R&    dE) :
        status(s),
        process(p), primary(is_p), particle(t), vtx0(v0), vtx1(v1), dE(dE), scorer_column(0),
        local_dE(0) {
        ;
    }

    /**
     * @brief Copy constructor.
     */
    CUDA_HOST_DEVICE
    track_t(const track_t& rhs) {
        scorer_column = rhs.scorer_column;
        status        = rhs.status;
        process       = rhs.process;
        particle      = rhs.particle;
        primary       = rhs.primary;
        vtx0          = rhs.vtx0;
        vtx1          = rhs.vtx1;
        dE            = rhs.dE;
        local_dE      = rhs.local_dE;
        c_node        = rhs.c_node;
        its           = rhs.its;
        ref_vector    = rhs.ref_vector;
    }

    /**
     * @brief Destructor.
     */
    CUDA_HOST_DEVICE
    ~track_t() {
        ;
    }

    /**
     * @brief Adds to the total energy deposited in this step.
     * @param e The amount of energy to deposit.
     */
    CUDA_HOST_DEVICE
    void
    deposit(R e) {
        dE += e;
    }

    /**
     * @brief Adds to the locally deposited energy in this step.
     * @param e The amount of local energy to deposit.
     */
    CUDA_HOST_DEVICE
    void
    local_deposit(R e) {
        local_dE += e;
    }

    /**
     * @brief Checks if the track's status is STOPPED.
     * @return `true` if the track is stopped, `false` otherwise.
     */
    CUDA_HOST_DEVICE
    bool
    is_stopped() {
        return status == STOPPED;
    }

    /**
     * @brief Reduces the length of the current step by a given ratio.
     * @param ratio The factor by which to shorten the step (0 < ratio < 1).
     */
    CUDA_HOST_DEVICE
    void
    shorten_step(R ratio)   //0 < ratio < 1
    {
        vtx1.pos = vtx0.pos + (vtx1.pos - vtx0.pos) * ratio;
    }

    /**
     * @brief Updates the post-step direction after a scattering event.
     * @param theta The polar scattering angle.
     * @param phi The azimuthal scattering angle.
     */
    CUDA_HOST_DEVICE
    void
    update_post_vertex_direction(const R& theta, const R& phi) {
        mqi::mat3x3<R> m_local(0, theta, phi);
        mqi::vec3<R>   d_local = m_local * ref_vector;   // rotate about the z-axis (dir)
        d_local.normalize();
        mqi::mat3x3<R> m_global(ref_vector, vtx1.dir);   // match dir to vtx1.dir
        vtx1.dir = m_global * d_local;
        vtx1.dir.normalize();
    }

    /**
     * @brief Updates the post-step position after a straight-line transport step.
     * @param len The length of the step.
     */
    CUDA_HOST_DEVICE
    void
    update_post_vertex_position(const R& len) {
        vtx1.pos = vtx0.pos + vtx0.dir * len;
    }

    /**
     * @brief Updates the post-step kinetic energy after an energy loss event.
     * @param e The amount of energy lost.
     */
    CUDA_HOST_DEVICE
    void
    update_post_vertex_energy(const R& e) {
        vtx1.ke -= e;
    }

    /**
     * @brief Finalizes the current step and prepares for the next one.
     * @details Sets the pre-step vertex (`vtx0`) to the current post-step vertex (`vtx1`) and resets energy deposition counters.
     */
    CUDA_HOST_DEVICE
    void
    move() {
        vtx0     = vtx1;
        dE       = 0;
        local_dE = 0;
    }

    /**
     * @brief Sets the track's status to STOPPED.
     */
    CUDA_HOST_DEVICE
    void
    stop() {
        status = STOPPED;
    }
};

/**
 * @brief A debugging function to assert the validity of a track's direction vectors.
 * @details Checks if the direction vectors in `vtx0` or `vtx1` contain NaN values or are zero vectors,
 * which would indicate a numerical error. If an error is found, it prints track details and exits.
 * @tparam R The floating-point type of the track.
 * @param trk The track to check.
 * @param id An optional identifier to print upon failure.
 */
template<typename R>
CUDA_HOST_DEVICE void
assert_track(const mqi::track_t<R>& trk, int8_t id = -1) {
    if (mqi::mqi_isnan(trk.vtx1.dir.x) || mqi::mqi_isnan(trk.vtx1.dir.y) ||
        mqi::mqi_isnan(trk.vtx1.dir.z)) {
        printf("id: %d\n", id);
        printf("vtx0.pos ");
        trk.vtx0.pos.dump();
        printf("vtx0.dir ");
        trk.vtx0.dir.dump();
        printf("vtx1.pos ");
        trk.vtx1.pos.dump();
        printf("vtx1.dir ");
        trk.vtx1.dir.dump();
        printf("There is nan in track direction\n");
        exit(1);
    } else if (trk.vtx0.dir.norm() < mqi::geometry_tolerance ||
               trk.vtx1.dir.norm() < mqi::geometry_tolerance) {
        printf("id: %d\n", id);
        printf("vtx0.pos ");
        trk.vtx0.pos.dump();
        printf("vtx0.dir ");
        trk.vtx0.dir.dump();
        printf("vtx1.pos ");
        trk.vtx1.pos.dump();
        printf("vtx1.dir ");
        trk.vtx1.dir.dump();
        printf("There is all-zeros in track direction\n");
        exit(1);
    }
}

}   // namespace mqi
#endif
