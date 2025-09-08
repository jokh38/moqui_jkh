#ifndef MQI_ROI_HPP
#define MQI_ROI_HPP
#include <moqui/base/mqi_common.hpp>

namespace mqi
{
/**
 * @enum roi_mapping_t
 * @brief Defines the method for mapping a transport grid index to a Region of Interest (ROI) scoring index.
 */
typedef enum
{
    DIRECT   = 0, ///< Direct mapping: The scoring index is the same as the transport grid index. Used when scoring the entire grid.
    INDIRECT = 1, ///< Indirect mapping: The scoring index is found through a lookup table (`start_[transport_index]`).
    CONTOUR  = 2  ///< Contour-based mapping: The scoring index is determined by searching through a run-length encoded contour map.
} roi_mapping_t;

/**
 * @class roi_t
 * @brief Manages the mapping from a global transport grid to a sparse Region of Interest (ROI) for scoring.
 * @details This class provides the logic to determine if a given voxel (identified by its transport grid index)
 * falls within an ROI and to find its corresponding index within the ROI's own flattened data array. It supports
 * different mapping strategies to handle various ROI definitions efficiently, especially on the GPU.
 */
class roi_t
{
public:
    roi_mapping_t method_;          ///< The mapping method to use (DIRECT, INDIRECT, or CONTOUR).
    uint32_t      original_length_; ///< The total number of voxels in the original, full transport grid.
    uint32_t      length_;          ///< The number of elements in the mapping arrays (`start_`, `stride_`, `acc_stride_`).

    // These array pointers are set externally, typically pointing to memory allocated on the GPU.
    uint32_t* start_;      ///< For CONTOUR: stores the starting transport index of each continuous segment (run). For INDIRECT: stores the mapped ROI index for each transport index.
    uint32_t* stride_;     ///< For CONTOUR: stores the length (number of consecutive voxels) of each segment.
    uint32_t* acc_stride_; ///< For CONTOUR: stores the accumulated length of segments, used to calculate the final ROI index.

public:
    /**
     * @brief Constructs an roi_t object.
     * @param m The mapping method.
     * @param n The total size of the original transport grid.
     * @param l The length of the mapping arrays.
     * @param s Pointer to the start/indirect-index array.
     * @param t Pointer to the stride array (for CONTOUR mapping).
     * @param a Pointer to the accumulated stride array (for CONTOUR mapping).
     */
    CUDA_HOST_DEVICE
    roi_t(roi_mapping_t m,
          uint32_t      n,
          int32_t       l = 0,
          uint32_t*     s = nullptr,
          uint32_t*     t = nullptr,
          uint32_t*     a = nullptr) :
        method_(m),
        original_length_(n), length_(l), start_(s), stride_(t), acc_stride_(a) {
        ;
    }

    /**
     * @brief Determines if a transport index `v` is inside the ROI.
     * @param v The transport grid index.
     * @return For CONTOUR, returns 1 if inside, -1 if outside. For others, returns the mapped index.
     */
    CUDA_HOST_DEVICE
    int32_t
    idx(const uint32_t& v) const {

        switch (method_) {
        case INDIRECT:
            return start_[v];
        case CONTOUR:
            return idx_contour(v);
        default:
            return v;
        }
    }

    /**
     * @brief Gets the final index within the flattened ROI mask for a given transport index.
     * @param v The transport grid index.
     * @return The corresponding index in the ROI's own data array, or -1 if outside the ROI.
     */
    CUDA_HOST_DEVICE
    int32_t
    get_mask_idx(const uint32_t& v) const {
        switch (method_) {
        case INDIRECT:
            return start_[v];
        case CONTOUR:
            return get_contour_idx(v);
        default:
            return v;
        }
    }

    /**
     * @brief Gets the total number of voxels in the ROI mask.
     * @return The size of the ROI.
     */
    CUDA_HOST_DEVICE
    int32_t
    get_mask_size() const {
        switch (method_) {
        case INDIRECT:
            return length_;
        case CONTOUR:
            return acc_stride_[length_ - 1];
        default:
            return original_length_;
        }
    }

    /**
     * @brief Calculates the ROI mask index for a transport index using the CONTOUR method.
     * @param v The transport grid index.
     * @return The final ROI index if `v` is within a contour segment, otherwise -1.
     */
    CUDA_HOST_DEVICE
    int32_t
    get_contour_idx(const uint32_t& v) const {
        int32_t  c        = this->lower_bound_cpp(v) - 1;
        uint32_t distance = v - start_[c];
        if (distance < stride_[c]) {
            /// is in stride
            if (c >= 1) distance += acc_stride_[c - 1];
            return distance;
        }
        return -1;   //invalid
    }

    /**
     * @brief Checks if a transport index is within any contour segment.
     * @param v The transport grid index.
     * @return 1 if `v` is inside a contour segment, otherwise -1.
     */
    CUDA_HOST_DEVICE
    int32_t
    idx_contour(const uint32_t& v) const {
        int32_t  c        = this->lower_bound_cpp(v) - 1;
        uint32_t distance = v - start_[c];
        if (distance < stride_[c]) {
            /// is in stride
            return 1;
        }
        return -1;   //invalid
    }

    /**
     * @brief Gets the mapped index for the DIRECT method (for completeness, though typically handled by default).
     * @param v The transport grid index.
     * @return The mapped ROI index from the lookup table.
     */
    CUDA_HOST_DEVICE
    int32_t
    idx_direct(const uint32_t& v) const {
        return start_[v];
    }

    /**
     * @brief A custom binary search implementation to find the lower bound of a value in the `start_` array.
     * @details This is a key performance component for the CONTOUR method, as it efficiently finds the
     * run-length-encoded segment that might contain the transport index `v`.
     * @param value The transport index to search for.
     * @return The index of the first element in `start_` that is greater than `value`.
     */
    CUDA_HOST_DEVICE
    int32_t
    lower_bound_cpp(const int32_t& value) const {
        int32_t first = 0;
        int32_t count = length_;
        int32_t step;

        int32_t it;
        while (count > 0) {

            it   = first;
            step = count / 2;
            it += step;

            if (start_[it] <= value) {
                first = ++it;
                count -= step + 1;
            } else
                count = step;
        }
        return first;
    }
};

}   // namespace mqi

#endif
