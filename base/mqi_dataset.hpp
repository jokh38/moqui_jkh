/// \file
///
/// \brief Defines a simplified interface for accessing DICOM datasets.
///
/// This file contains the definition of the `mqi::dataset` class, which provides
/// a wrapper around the `gdcm::DataSet` to simplify accessing and retrieving
/// values from DICOM files, especially those with nested sequences.

#ifndef MQI_DATASET_H
#define MQI_DATASET_H

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "gdcmByteValue.h"
#include "gdcmDataElement.h"
#include "gdcmDataSet.h"
#include "gdcmDicts.h"
#include "gdcmGlobal.h"
#include "gdcmItem.h"
#include "gdcmSequenceOfItems.h"
#include "gdcmVM.h"
#include "gdcmVR.h"

namespace mqi
{

/// \enum modality_type
/// \brief Enumerates DICOM modality types based on the Modality (0x0008, 0x0060) tag.
typedef enum
{
    RTPLAN,      ///< Radiation Therapy Plan
    IONPLAN,     ///< Ion Therapy Plan
    RTRECORD,    ///< Radiation Therapy Treatment Record
    IONRECORD,   ///< Ion Therapy Treatment Record
    RTIMAGE,     ///< Radiation Therapy Image (e.g., CT)
    RTSTRUCT,    ///< Radiation Therapy Structure Set
    RTDOSE,      ///< Radiation Therapy Dose
    UNKNOWN_MOD  ///< Unknown or unspecified modality
} modality_type;

/// \struct beam_id_type
/// \brief Represents a beam identifier, which can be either a number or a string.
struct beam_id_type {
    /// \enum
    /// \brief The type of the beam identifier.
    enum
    {
        NUM, ///< The identifier is a number.
        STR  ///< The identifier is a string.
    } type;
    /// \union
    /// \brief The value of the beam identifier.
    union {
        int         number; ///< The numeric identifier.
        const char* name;   ///< The string identifier.
    };
};

/// \brief A map defining sequence tags for different DICOM modalities.
///
/// This map provides a convenient way to look up the `gdcm::Tag` for common
/// sequences (like beam or control point sequences) based on the modality
/// of the DICOM file (e.g., RTPLAN, IONPLAN).
static const std::map<const modality_type, const std::map<const std::string, const gdcm::Tag>>
  seqtags_per_modality = {
      //modality, seq_name, tag
      { RTPLAN,
        {
          { "beam", gdcm::Tag(0x300a, 0x00b0) }
          //{"wedge", gdcm::Tag()},
          //{"mlc"  , gdcm::Tag()},
        } },
      {
        IONPLAN,
        { { "beam", gdcm::Tag(0x300a, 0x03a2) },
          { "snout", gdcm::Tag(0x300a, 0x030c) },
          { "rs", gdcm::Tag(0x300a, 0x0314) },     ///< range shifter sequence
          { "rsss", gdcm::Tag(0x300a, 0x0360) },   ///< range shifter setting sequence
          { "blk", gdcm::Tag(0x300a, 0x03a6) },
          { "comp", gdcm::Tag(0x300a, 0x0e2a) },
          { "ctrl", gdcm::Tag(0x300a, 0x03a8) } }   ///< ion control point sequence
      },
      { IONRECORD,
        {
          { "beam", gdcm::Tag(0x3008, 0x0021) },
          { "snout", gdcm::Tag(0x3008, 0x00f0) },
          { "rs", gdcm::Tag(0x3008, 0x00f2) },
          { "rsss", gdcm::Tag(0x300a, 0x0360) },   ///< range shifter setting sequence
          { "blk", gdcm::Tag(0x3008, 0x00d0) },
          { "comp", gdcm::Tag(0x3008, 0x00c0) },
          { "ctrl", gdcm::Tag(0x3008, 0x0041) },
          { "machine", gdcm::Tag(0x300a, 0x0206) }   ///< treatment machine record
        } }
  };

/// \class dataset
/// \brief A wrapper for `gdcm::DataSet` to simplify access to DICOM data elements and sequences.
///
/// The `dataset` class provides an intuitive interface for querying DICOM data.
/// It allows accessing data elements by tag or keyword using the `[]` operator,
/// and accessing nested datasets within sequences using the `()` operator.
/// It also includes methods for converting DICOM data values into standard C++ types.
///
/// Example Usage:
/// ```cpp
/// // Get a value by keyword
/// std::vector<int> num_blocks;
/// block_ds->get_values("NumberOfBlocks", num_blocks);
///
/// // Access a sequence by keyword
/// const std::vector<const dataset*> beam_datasets = plan_ds("IonBeamSequence");
/// ```
class dataset
{
protected:
    ///< A lookup table for nested datasets, mapping tag and keyword to a vector of child datasets.
    std::vector<std::tuple<const gdcm::Tag, const std::string, std::vector<const dataset*>>>
      ds_lut_;

    ///< The underlying GDCM dataset.
    const gdcm::DataSet gdcm_ds_;

public:
    /// \brief Constructs a dataset object from a `gdcm::DataSet`.
    /// \param[in] d The `gdcm::DataSet` to wrap.
    /// \param[in] include_sq If true, recursively processes sequences to build the `ds_lut_`.
    dataset(const gdcm::DataSet& d, const bool include_sq = true) : gdcm_ds_(d) {
        const gdcm::Dicts& dicts = gdcm::Global::GetInstance().GetDicts();

        if (!include_sq) { return; }

        for (auto el = gdcm_ds_.Begin(); el != gdcm_ds_.End(); ++el) {
            /// el->getValue() doesn't guarantee it has value.
            const gdcm::Tag&       tag   = el->GetTag();
            const gdcm::DictEntry& entry = dicts.GetDictEntry(tag);

            if (!(entry.GetVR() & gdcm::VR::SQ)) continue;

            gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = el->GetValueAsSQ();
            std::vector<const dataset*>               tmp(0);
            for (size_t i = 1; i <= sqi->GetNumberOfItems(); ++i) {
                const gdcm::Item& itm = sqi->GetItem(i);
                tmp.push_back(new mqi::dataset(itm.GetNestedDataSet()));
            }
            ds_lut_.push_back(std::make_tuple(tag, entry.GetKeyword(), tmp));

        }   //gdcm_ds_
    }

    /// \brief Accesses a `gdcm::DataElement` by its tag.
    /// \param[in] t The `gdcm::Tag` of the data element.
    /// \return A constant reference to the `gdcm::DataElement`. Returns an empty element if not found.
    const gdcm::DataElement&
    operator[](const gdcm::Tag& t) const {
        if (gdcm_ds_.FindDataElement(t) && !gdcm_ds_.GetDataElement(t).IsEmpty()) {
            return gdcm_ds_.GetDataElement(t);
        } else {
            static gdcm::DataElement empty;
            empty.Clear();   //invalid VR & VL
            return empty;
        }
    }

    /// \brief Accesses a `gdcm::DataElement` by its keyword.
    /// \param[in] keyword The keyword of the data element as a C-style string.
    /// \return A constant reference to the `gdcm::DataElement`. Returns an empty element if not found.
    const gdcm::DataElement&
    operator[](const char* keyword) const {
        std::string        tmp(keyword);
        const gdcm::Dicts& dicts = gdcm::Global::GetInstance().GetDicts();

        for (auto it = gdcm_ds_.Begin(); it != gdcm_ds_.End(); ++it) {
            const gdcm::DictEntry& entry = dicts.GetDictEntry(it->GetTag());
            if (!tmp.compare(entry.GetKeyword())) return (*it);
        }
        static gdcm::DataElement empty;
        empty.Clear();   //invalid VR & VL
        return empty;
    }

    /// \brief Accesses a sequence of nested datasets by its tag.
    /// \param[in] t The `gdcm::Tag` of the sequence.
    /// \return A vector of pointers to the constant `dataset` objects within the sequence. Returns an empty vector if not found.
    std::vector<const dataset*>
    operator()(const gdcm::Tag& t) const {
        for (auto& i : ds_lut_)
            if (std::get<0>(i) == t) { return std::get<2>(i); }
        std::vector<const dataset*> empty(0);
        return empty;
    }

    /// \brief Accesses a sequence of nested datasets by its keyword.
    /// \param[in] s The keyword of the sequence as a C-style string.
    /// \return A vector of pointers to the constant `dataset` objects within the sequence. Returns an empty vector if not found.
    std::vector<const dataset*>
    operator()(const char* s) const {
        std::string tmp(s);
        for (auto& i : ds_lut_)
            if (!tmp.compare(std::get<1>(i))) return std::get<2>(i);

        std::vector<const dataset*> empty(0);
        return empty;
    }

    /// \brief Destructor.
    ///
    /// Cleans up the dynamically allocated nested `dataset` objects.
    ~dataset() {
        for (auto& i : ds_lut_) {
            for (auto& j : std::get<2>(i)) {
                delete j;
            }
            std::get<2>(i).clear();
        }
        ds_lut_.clear();
    }

    /// \brief Prints all data elements in the dataset to the console.
    ///
    /// This method is useful for debugging and inspecting the content of a dataset.
    void
    dump() const {
        const gdcm::Dicts& dicts = gdcm::Global::GetInstance().GetDicts();

        for (auto el = gdcm_ds_.Begin(); el != gdcm_ds_.End(); ++el) {
            const gdcm::DictEntry& entry = dicts.GetDictEntry(el->GetTag());
            std::cout << entry.GetKeyword() << ": " << el->GetTag() << ", VR: " << el->GetVR()
                      << ", length: " << el->GetVL() << std::endl;
        }

        std::cout << " LUT: --- size: " << ds_lut_.size() << std::endl;
        for (auto& i : ds_lut_) {
            std::cout << "     (" << std::get<0>(i) << ", " << std::get<1>(i) << ", "
                      << std::get<2>(i).size() << " )";
        }
        std::cout << std::endl;
    }

    /// \brief Extracts values from a `gdcm::ByteValue` into a vector of floats.
    /// \param[in] bv Pointer to the `gdcm::ByteValue` containing the data.
    /// \param[in] vr The Value Representation (VR) of the data.
    /// \param[in] vm The Value Multiplicity (VM) of the data.
    /// \param[in] vl The Value Length (VL) of the data.
    /// \param[out] res A reference to a vector of floats to store the results.
    void
    get_values(const gdcm::ByteValue* bv,
               const gdcm::VR         vr,
               const gdcm::VM         vm,
               const gdcm::VL         vl,
               std::vector<float>&    res) const {
        res.clear();

        if (vr & gdcm::VR::VRBINARY) {
            assert(vr & gdcm::VR::FL);

            size_t ndim = vl / sizeof(float);   ///< T should be given 'float'
            res.resize(ndim);
            bv->GetBuffer((char*) &res[0], ndim * sizeof(float));
        } else if (vr & gdcm::VR::VRASCII) {
            //ascii & numeric (int, float)
            //ascii & IS, DS, -> decimal strings...
            std::string       s    = std::string(bv->GetPointer(), bv->GetLength());
            size_t            beg  = 0;
            size_t            next = std::string::npos;
            const std::string tok("\\");
            do {
                next = s.find_first_of(tok, beg);
                res.push_back(std::stof(s.substr(beg, next)));
                beg = next + tok.size();
            } while (next != std::string::npos);
        }
    }

    /// \brief Extracts values from a `gdcm::ByteValue` into a vector of integers.
    /// \param[in] bv Pointer to the `gdcm::ByteValue` containing the data.
    /// \param[in] vr The Value Representation (VR) of the data.
    /// \param[in] vm The Value Multiplicity (VM) of the data.
    /// \param[in] vl The Value Length (VL) of the data.
    /// \param[out] res A reference to a vector of integers to store the results.
    /// \param[in] show (Unused) A boolean flag.
    void
    get_values(const gdcm::ByteValue* bv,
               const gdcm::VR         vr,
               const gdcm::VM         vm,
               const gdcm::VL         vl,
               std::vector<int>&      res,
               bool                   show = false) const {
        //works for VR of FL, DS, ...
        res.clear();

        size_t ndim = vm.GetLength();

        if (ndim == 0) {
            if (vr & gdcm::VR::FL) {       //
                ndim = vl / sizeof(int);   //T should be given 'float'
            }
        }
        assert(ndim >= 1);   //From here, dim shouldn't be 0.a

        if (vr & gdcm::VR::VRBINARY) {
            res.resize(ndim);
            bv->GetBuffer((char*) &res[0], ndim * sizeof(int));
        } else if (vr & gdcm::VR::VRASCII) {
            std::string s = std::string(bv->GetPointer(), bv->GetLength());

            size_t            beg  = 0;
            size_t            next = std::string::npos;
            const std::string tok("\\");
            do {
                next = s.find_first_of(tok, beg);
                res.push_back(std::stoi(s.substr(beg, next)));
                beg = next + tok.size();
            } while (next != std::string::npos);
        }
    }

    /// \brief Extracts values from a `gdcm::ByteValue` into a vector of strings.
    /// \param[in] bv Pointer to the `gdcm::ByteValue` containing the data.
    /// \param[in] vr The Value Representation (VR) of the data.
    /// \param[in] vm The Value Multiplicity (VM) of the data.
    /// \param[in] vl The Value Length (VL) of the data.
    /// \param[out] res A reference to a vector of strings to store the results.
    void
    get_values(const gdcm::ByteValue*    bv,
               const gdcm::VR            vr,
               const gdcm::VM            vm,
               const gdcm::VL            vl,
               std::vector<std::string>& res) const {
        //works for VR of FL, DS, ...
        res.clear();

        size_t ndim = vm.GetLength();

        assert(ndim >= 1);   //From here, dim shouldn't be 0.

        if (vr & gdcm::VR::VRASCII) {
            std::string s = std::string(bv->GetPointer(), bv->GetLength());

            size_t            beg  = 0;
            size_t            next = std::string::npos;
            const std::string tok("\\");
            do {
                next = s.find_first_of(tok, beg);
                res.push_back(s.substr(beg, next));
                beg = next + tok.size();
            } while (next != std::string::npos);
        }
    }

    /// \brief A template method to get values for a given keyword.
    /// \tparam T The type of the values to retrieve (e.g., int, float, std::string).
    /// \param[in] keyword The keyword of the data element.
    /// \param[out] result A reference to a vector of type T to store the results.
    template<typename T>
    void
    get_values(const char* keyword, std::vector<T>& result) const {
        result.clear();
        const gdcm::DataElement& el = (*this)[keyword];
        if (el.IsEmpty()) {
            result.resize(0);
            return;
        }
        const gdcm::Dicts&     dicts = gdcm::Global::GetInstance().GetDicts();
        const gdcm::Tag&       tag   = el.GetTag();
        const gdcm::DictEntry& entry = dicts.GetDictEntry(tag);

        this->get_values(el.GetByteValue(), entry.GetVR(), entry.GetVM(), el.GetVL(), result);
    }

    /// \brief A template method to get values from a `gdcm::DataElement`.
    /// \tparam T The type of the values to retrieve.
    /// \param[in] el The `gdcm::DataElement` from which to extract values.
    /// \param[out] result A reference to a vector of type T to store the results.
    template<typename T>
    void
    get_values(const gdcm::DataElement& el, std::vector<T>& result) const {
        result.clear();
        if (el.IsEmpty()) {
            result.resize(0);
            return;
        }
        const gdcm::Dicts&     dicts = gdcm::Global::GetInstance().GetDicts();
        const gdcm::Tag&       tag   = el.GetTag();
        const gdcm::DictEntry& entry = dicts.GetDictEntry(tag);

        this->get_values(el.GetByteValue(), entry.GetVR(), entry.GetVM(), el.GetVL(), result);
    }

    /// \brief Prints information about a `gdcm::DataElement`.
    /// \param[in] el The data element to print.
    ///
    /// This method is useful for debugging, showing the keyword, VM, VR, VL, and value of a data element.
    void
    print_dataelement(const gdcm::DataElement& el) const {
        const gdcm::Dicts& dicts = gdcm::Global::GetInstance().GetDicts();

        //convert ByteValue by using GetVR/VM
        const gdcm::ByteValue* bv = el.GetByteValue();

        const gdcm::Tag        tag   = el.GetTag();
        const gdcm::DictEntry& entry = dicts.GetDictEntry(tag);
        //el.GetType()
        //VM1-n doesn't have a number
        std::cout << "---- " << entry.GetKeyword() << " VM:" << entry.GetVM().GetLength()
                  << " VR:" << entry.GetVR() << ", GetVL: " << el.GetVL()
                  << ", GetByteValues: " << el.GetByteValue() << ", VR length"
                  << el.GetVR().GetLength() << ", VR size"
                  << el.GetVR().GetVRString(
                       el.GetVR())   //don't call (.GetSize of VR) this. it produced error
                  << ",    value:";
        if (el.IsEmpty()) {
            std::cout << "no value." << std::endl;
        } else {
            if (entry.GetVR() & gdcm::VR::FL) {

                const size_t ndim = el.GetVL() / sizeof(float);

                float c[ndim];
                bv->GetBuffer((char*) c, sizeof(c));

                std::cout << "FL encoding: ";
                for (size_t i = 0; i < ndim; ++i) {
                    std::cout << "   " << i << ": " << c[i];
                }
                std::cout << std::endl;

            } else {
            }
        }
    }
};

}   // namespace mqi
#endif
