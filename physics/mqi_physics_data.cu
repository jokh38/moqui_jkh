#include "mqi_physics_data.hpp"
#include <iostream>
#include <stdexcept>
#include <cstring>

#include <moqui/base/mqi_p_ionization.hpp>
#include <moqui/base/mqi_pp_elastic.hpp>
#include <moqui/base/mqi_po_elastic.hpp>
#include <moqui/base/mqi_po_inelastic.hpp>

namespace mqi {

physics_data_manager::physics_data_manager() {
    // Constructor
}

physics_data_manager::~physics_data_manager() {
    // Cleanup: destroy texture objects and free CUDA arrays
    for (auto const& [key, val] : texture_objects_) {
        cudaDestroyTextureObject(val);
    }
    for (auto const& [key, val] : cuda_arrays_) {
        cudaFreeArray(val);
    }
}

void physics_data_manager::create_texture(const std::string& data_type, int width, int height, const float* h_data) {
    // Create channel description
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

    // Allocate CUDA array
    cudaArray* cu_array;
    cudaMallocArray(&cu_array, &channelDesc, width, height);
    cuda_arrays_[data_type] = cu_array;

    // Copy data to CUDA array
    cudaMemcpy2DToArray(cu_array, 0, 0, h_data, width * sizeof(float), width * sizeof(float), height, cudaMemcpyHostToDevice);

    // Create texture object
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = cu_array;

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.addressMode[1] = cudaAddressModeClamp;
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 0;

    cudaTextureObject_t texObject = 0;
    cudaCreateTextureObject(&texObject, &resDesc, &texDesc, NULL);
    texture_objects_[data_type] = texObject;
}


void physics_data_manager::initialize() {
    // For now, let's assume the tables are 1D and have 600 elements.
    // The 2D texture will have a height of 1.
    // I am creating one large texture for all cross sections.
    // The y-coordinate will be used to select the cross section type.

    const int width = 600;
    const int height = 4; // for p_ion, pp_e, po_e, po_i

    std::vector<float> h_cross_sections(width * height);
    memcpy(h_cross_sections.data() + width * 0, cs_p_ion_table, width * sizeof(float));
    memcpy(h_cross_sections.data() + width * 1, cs_pp_e_g4_table, width * sizeof(float));
    memcpy(h_cross_sections.data() + width * 2, cs_po_e_g4_table, width * sizeof(float));
    memcpy(h_cross_sections.data() + width * 3, cs_po_i_g4_table, width * sizeof(float));

    create_texture("cross_section", width, height, h_cross_sections.data());

    // Stopping power is a separate texture
    create_texture("stopping_power", width, 1, restricted_stopping_power_table);
}

cudaTextureObject_t physics_data_manager::get_texture_object(const std::string& data_type) const {
    auto it = texture_objects_.find(data_type);
    if (it == texture_objects_.end()) {
        throw std::runtime_error("Texture object not found for data type: " + data_type);
    }
    return it->second;
}

} // namespace mqi
