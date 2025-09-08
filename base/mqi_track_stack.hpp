
#ifndef MQI_TRACK_STACK_HPP
#define MQI_TRACK_STACK_HPP

#include <moqui/base/mqi_node.hpp>
#include <moqui/base/mqi_track.hpp>

namespace mqi
{

/**
 * @class track_stack_t
 * @brief A simple, fixed-size stack for managing secondary particle tracks.
 * @details This class provides a Last-In, First-Out (LIFO) container for `track_t` objects.
 * In a Monte Carlo simulation, when a particle interaction creates new secondary particles,
 * they are pushed onto this stack. After the primary particle's transport is complete,
 * tracks are popped from the stack to be transported. This ensures all particle lineages are
 * fully simulated. The stack has a fixed capacity, which is larger in debug builds to
 * facilitate more detailed tracking.
 * @tparam R The floating-point type for the track data (e.g., float, double).
 */
template<typename R>
class track_stack_t
{

public:
#ifdef __PHYSICS_DEBUG__
    const uint16_t limit = 200;   ///< The maximum number of tracks the stack can hold in debug mode.
    track_t<R>     tracks[200]; ///< The underlying array to store track objects in debug mode.
#else
    const uint16_t limit = 10;    ///< The maximum number of tracks the stack can hold in release mode.
    track_t<R>     tracks[10];  ///< The underlying array to store track objects in release mode.
#endif
    uint16_t idx = 0; ///< The current number of tracks in the stack, acting as a stack pointer. An index of 0 means the stack is empty.

    /**
     * @brief Default constructor.
     */
    CUDA_HOST_DEVICE
    track_stack_t() {
        ;
    }

    /**
     * @brief Destructor.
     */
    CUDA_HOST_DEVICE
    ~track_stack_t() {
        ;
    }

    /**
     * @brief Pushes a secondary track onto the top of the stack.
     * @details If the stack is already full (i.e., `idx >= limit`), the operation is ignored.
     * @param trk The track to be added to the stack.
     */
    CUDA_HOST_DEVICE
    void
    push_secondary(const track_t<R>& trk) {

        if (idx < limit) {
            tracks[idx] = trk;
            ++idx;
        }
    }

    /**
     * @brief Checks if the stack is empty.
     * @return `true` if the stack contains no tracks, `false` otherwise.
     */
    CUDA_HOST_DEVICE
    bool
    is_empty(void) {
        return idx == 0;
    }

    /**
     * @brief Removes and returns the track from the top of the stack.
     * @return A copy of the track that was at the top of the stack.
     */
    CUDA_HOST_DEVICE
    track_t<R>
    pop(void) {
        ///copy
        return tracks[--idx];
    }

    /**
     * @brief Provides direct array-like access to the tracks in the stack.
     * @param i The index of the track to access.
     * @return A reference to the track at the specified index.
     */
    CUDA_HOST_DEVICE
    track_t<R>&
    operator[](uint16_t i) {
        return tracks[i];
    }
};

}   // namespace mqi
#endif
