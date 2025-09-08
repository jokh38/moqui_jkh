/**
 * @file
 * @brief Defines the node structure for the geometry hierarchy.
 */
#ifndef MQI_NODE_HPP
#define MQI_NODE_HPP

#include <moqui/base/mqi_grid3d.hpp>
#include <moqui/base/mqi_scorer.hpp>

#if defined(__CUDACC__)
#include <cuda_fp16.h>
#endif

namespace mqi
{

/**
 * @struct node_t
 * @brief Represents a node in a hierarchical geometry, containing geometry, scorers, and child nodes.
 * @tparam R The floating-point type for scoring and geometry values (e.g., float or double).
 * @details This structure is a core component of the simulation geometry, allowing for nested or complex arrangements. Each node can have its own geometry definition and a set of scorers to record data.
 */
template<typename R>
struct node_t {
    grid3d<mqi::density_t, R>* geo = nullptr;   ///< A pointer to the node's geometry, defined as a 3D grid.

    uint16_t         n_scorers = 0;         ///< The number of scorers attached to this node.
    scorer<R>**      scorers   = nullptr;   ///< An array of pointers to the scorer objects.
    mqi::key_value** scorers_data =
        nullptr;   ///< An array of pointers to the primary data for each scorer (e.g., dose).
    mqi::key_value** scorers_count =
        nullptr;   ///< An array of pointers to the count data for each scorer (for statistical calculations).
    mqi::key_value** scorers_mean =
        nullptr;   ///< An array of pointers to the mean value data for each scorer.
    mqi::key_value** scorers_variance =
        nullptr;   ///< An array of pointers to the variance data for each scorer.

    uint16_t           n_children = 0;      ///< The number of child nodes.
    struct node_t<R>** children   = nullptr;   ///< An array of pointers to the child nodes.
};

}   // namespace mqi
#endif
