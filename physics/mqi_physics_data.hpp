#ifndef MQI_PHYSICS_DATA_HPP
#define MQI_PHYSICS_DATA_HPP

#include <cuda_runtime.h>
#include <vector>
#include <string>
#include <map>

namespace mqi {

///< A class to manage physics data for GPU transport kernels.
// It loads data from files, creates CUDA arrays, and binds them to texture objects.
class physics_data_manager {
public:
    ///< Default constructor
    physics_data_manager();

    ///< Destructor
    ~physics_data_manager();

    ///< Load data and initialize textures
    void initialize();

    ///< Get the texture object for a given data type (e.g., "stopping_power")
    cudaTextureObject_t get_texture_object(const std::string& data_type) const;

private:
    ///< Creates a 2D CUDA array and texture object
    void create_texture(const std::string& data_type, int width, int height, const float* h_data);

    ///< A map to store texture objects, keyed by data type
    std::map<std::string, cudaTextureObject_t> texture_objects_;

    ///< A map to store cuda arrays, keyed by data type
    std::map<std::string, cudaArray*> cuda_arrays_;
};

} // namespace mqi

#endif // MQI_PHYSICS_DATA_HPP
