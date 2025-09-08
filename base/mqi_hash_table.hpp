#ifndef MQI_HASH_TABLE_CLASS_HPP
#define MQI_HASH_TABLE_CLASS_HPP

#include <cstring>
#include <moqui/base/mqi_common.hpp>

namespace mqi
{

/**
 * @struct key_value
 * @brief Represents a key-value pair for a hash table.
 * @details This structure holds two keys and a corresponding double-precision value. It is used as the basic element in the hash table implementation.
 */
struct key_value {
    mqi::key_t key1;  ///< The first key component.
    mqi::key_t key2;  ///< The second key component.
    double     value; ///< The value associated with the key pair.
};

/**
 * @brief Initializes a hash table on the host.
 * @param[out] table A pointer to the array of key_value structs to be initialized.
 * @param[in] max_capacity The maximum capacity of the hash table.
 * @details This function initializes all entries in the table, setting keys to `mqi::empty_pair` and values to 0.
 */
void
init_table(key_value* table, uint32_t max_capacity) {
    //// Multithreading?
    for (int i = 0; i < max_capacity; i++) {
        table[i].key1  = mqi::empty_pair;
        table[i].key2  = mqi::empty_pair;
        table[i].value = 0;
    }
}

/**
 * @brief Initializes a hash table on a CUDA device.
 * @tparam R The floating-point type (e.g., float or double).
 * @param[out] table A pointer to the device memory where the key_value structs are stored.
 * @param[in] max_capacity The maximum capacity of the hash table.
 * @details This CUDA kernel initializes the `value` field of each entry in the table to 0. Note: This function only initializes the value, not the keys.
 */
template<typename R>
CUDA_GLOBAL void
init_table_cuda(key_value* table, uint32_t max_capacity) {

    //// Multithreading?
    for (int i = 0; i < max_capacity; i++) {
        table[i].value = 0;
    }
    //#endif
}

/**
 * @brief A test function to print a specific entry from the hash table on a CUDA device.
 * @tparam R The floating-point type (e.g., float or double).
 * @param[in] data A pointer to the device memory where the key_value structs are stored.
 * @details This CUDA kernel prints the key and value of a hardcoded index for debugging purposes.
 */
template<typename R>
CUDA_GLOBAL void
test_print(mqi::key_value* data) {
    uint32_t ind = 512 * 512 * 200 * 4 - 1;
    printf("data[0].key1 %d data[0].key2 %d data[0].value %d\n",
           data[ind].key1,
           data[ind].key2,
           data[ind].value);
}
}   // namespace mqi
#endif
