#ifndef MQI_TRANSPORT_EVENT_HPP
#define MQI_TRANSPORT_EVENT_HPP

#include <moqui/base/mqi_track.hpp>
#include <moqui/base/mqi_physics_constants.hpp>
#include <moqui/base/mqi_vec.hpp>
#include <moqui/base/mqi_grid3d.hpp>

namespace mqi {

///
/// @brief Event-by-event particle transport kernel using Woodcock tracking.
///
/// This kernel transports particles one interaction at a time, offering higher
/// accuracy at the cost of computational intensity compared to condensed history.
/// It uses Woodcock tracking to avoid costly boundary crossings in heterogeneous media.
///
/// @param tracks           Device pointer to the array of particle tracks.
/// @param n_particles      The total number of particles to transport.
/// @param patient          Device pointer to the patient geometry grid.
/// @param xsec_texture     The texture object for accessing cross-section data.
/// @param stop_pow_texture The texture object for accessing stopping power data.
/// @param max_sigma        The maximum total cross-section used for Woodcock tracking.
///
template <typename T>
__global__ void
transport_event_by_event_kernel(mqi::thrd_t*      d_threads,
                                mqi::node_t<T>*   world,
                                mqi::vertex_t<T>* vertices,
                                int               n_particles,
                                uint32_t*         d_tracked_particles,
                                cudaTextureObject_t stop_pow_texture,
                                cudaTextureObject_t xsec_texture,
                                float             max_sigma)
{
    // Thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_particles) return;

    // Shared memory stack for secondary particles
    __shared__ track_t<T> secondary_stack[256];
    __shared__ int stack_top;
    if (threadIdx.x == 0) {
        stack_top = -1;
    }
    __syncthreads();

    // Get the primary particle for this thread
    mqi::mqi_rng* rng = &d_threads[idx].rnd_generator;
    track_t<T> p(vertices[idx]);

    // Main transport loop for the particle
    while (!p.is_stopped()) {
        // 1. Woodcock Tracking: Determine distance to next potential interaction site
        // - Sample a random number 'r'
        // - distance = -mqi_ln(r) / max_sigma
        // - Move particle by 'distance'
        // - At the new position, get the material and its actual total cross-section (sigma_actual)
        // - Sample another random number 'r2'
        // - If r2 < (sigma_actual / max_sigma), a real interaction occurs.
        // - Otherwise, it's a "null-collision" and the particle continues straight.

        // Placeholder for Woodcock tracking logic
        // p.pos += p.dir * distance;


        // 2. If a real interaction occurs:
        // - Sample the interaction type (e.g., ionization, elastic scatter)
        // - Update particle energy, direction, etc.
        // - If a secondary particle is created:
        //      - Atomically increment stack_top
        //      - secondary_stack[stack_top] = new_secondary_particle;

        // Placeholder for interaction physics
        p.vtx0.ke -= 0.1; // Dummy energy loss


        // 3. Check for particle death (energy below cutoff, left geometry)
        if (p.vtx0.ke <= 0.1) {
            p.stop();
        }

        // 4. Process secondary particles from the stack
        // - If this particle dies, and the stack is not empty, pop a secondary
        //   and start transporting it.
        // - A sync mechanism is needed to ensure all threads in a block
        //   contribute to and pull from the same stack correctly.
        if (p.is_stopped()) {
            __syncthreads(); // Sync before accessing stack
            int current_top = atomicSub(&stack_top, 1);
            if (current_top >= 0) {
                p = secondary_stack[current_top];
            }
        }

        // Advanced Optimizations (placeholders):
        // - __shfl_sync() could be used for fast data exchange within a warp.
        // - Dynamic parallelism: if the secondary stack grows too large, a new
        //   child kernel could be launched to process it.
        // - Intelligent energy cutoff: cutoff value could be dependent on the material
        //   at the particle's current position.

    } // end while(!p.is_stopped())

    // Write final state of the particle back to global memory
    // tracks[idx] = p; // 'tracks' is not defined in the kernel signature
}

} // namespace mqi

#endif // MQI_TRANSPORT_EVENT_HPP
