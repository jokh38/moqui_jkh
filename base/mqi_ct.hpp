/// \file
///
/// \brief Defines a class for handling 3D CT image data.
///
/// This file provides the `ct` class, which represents a 3D image grid with
/// 16-bit integer pixels. It is designed to read and process CT data from
/// DICOM files, managing the grid structure, pixel data, and relevant metadata.

#ifndef MQI_CT_H
#define MQI_CT_H

#include <sys/stat.h>

#include "gdcmAttribute.h"
#include "gdcmDirectory.h"
#include "gdcmIPPSorter.h"
#include "gdcmImageReader.h"
#include "gdcmScanner.h"

#include <moqui/base/mqi_matrix.hpp>
#include <moqui/base/mqi_rect3d.hpp>

namespace mqi
{

/// \class ct
/// \brief Represents a 3D CT image grid with 16-bit integer pixels.
/// \tparam R The floating-point type for grid coordinates (e.g., float, double).
///
/// This class inherits from `rect3d` and is specialized for handling CT data.
/// It includes functionality to read a directory of DICOM files, sort them,
/// extract metadata, and load the pixel data into a 3D grid.
template<typename R>
class ct : public rect3d<int16_t, R>
{
protected:
    ///< A vector of DICOM image file paths, sorted in ascending Z-order.
    std::vector<std::string> files_;

    ///< A map from SOP Instance UID to the corresponding file path.
    std::map<std::string, std::string> uid2file_;

    ///< The directory containing the CT files.
    char* ct_dir;

    ///< The pixel spacing in the X-dimension.
    R dx_;
    ///< The pixel spacing in the Y-dimension.
    R dy_;
    ///< An array of pixel spacings in the Z-dimension to handle variable slice thickness.
    R* dz_;

public:
    /// \brief Default constructor.
    CUDA_HOST
    ct() {
        ;
    }

    /// \brief Constructs a ct object and initializes its structure from a DICOM directory.
    /// \param[in] f The path to the directory containing the CT DICOM files.
    /// \param[in] is_print If true, prints the detected file information to the console.
    /// \note This constructor sets up the grid dimensions and properties but does not
    ///       load the pixel data. Call `load_data()` to read the pixel values.
    CUDA_HOST
    ct(std::string f, bool is_print = false) {
        ct_dir = new char[f.length() + 1];
        strcpy(ct_dir, f.c_str());
        //http://gdcm.sourceforge.net/html/SortImage_8cxx-example.html#_a5
        gdcm::Directory dir;
        dir.Load(this->ct_dir, false);   ///< non-recursive

        const gdcm::Directory::FilenamesType& all_files = dir.GetFilenames();

        gdcm::Scanner   scanner;
        const gdcm::Tag ct_tag = gdcm::Tag(0x08, 0x60);
        scanner.AddTag(ct_tag);
        scanner.Scan(all_files);

        auto ct_files = scanner.GetAllFilenamesFromTagToValue(ct_tag, "CT");

        gdcm::IPPSorter ippsorter;
        ippsorter.SetComputeZSpacing(false);
        ippsorter.Sort(ct_files);   ///< asending along z
        files_ = ippsorter.GetFilenames();

        gdcm::Scanner s;
        s.AddTag(gdcm::Tag(0x0008, 0x0018));   ///< SOP instance UID
        s.AddTag(gdcm::Tag(0x0020, 0x0032));   ///< Image Position (Patient)
        s.AddTag(gdcm::Tag(0x0028, 0x0010));   ///< Rows
        s.AddTag(gdcm::Tag(0x0028, 0x0011));   ///< Columns
        s.AddTag(gdcm::Tag(0x0028, 0x0030));   ///< Pixel spacing
        s.AddTag(gdcm::Tag(0x0018, 0x0050));   ///< Slice thickness

        if (!s.Scan(files_)) assert("scan fail.");

        size_t nx;                   ///< columns
        size_t ny;                   ///< rows
        size_t nz = files_.size();   ///< number of images.

        rect3d<int16_t, R>::z_     = new R[nz];
        rect3d<int16_t, R>::dim_.z = nz;
        double x0, y0;
        dz_ = new R[nz];
        float dz;
        for (size_t i = 0; i < nz; ++i) {
            gdcm::Scanner::TagToValue const& m0 = s.GetMapping(files_[i].c_str());

            std::string  img_position(m0.find(gdcm::Tag(0x0020, 0x0032))->second);
            unsigned int deli0        = img_position.find_first_of('\\');
            unsigned int deli1        = img_position.find_last_of('\\');
            rect3d<int16_t, R>::z_[i] = (R) (std::stod(img_position.substr(deli1 + 1)));
            dz_[i]                    = std::stod(m0.find(gdcm::Tag(0x0018, 0x0050))->second);

            ///< We only determine rows, colums, x0, y0, dx, and dy with first image
            if (i == 0) {
                x0 = std::stod(img_position.substr(0, deli0));
                y0 = std::stod(img_position.substr(deli0 + 1, deli1 - deli0));

                ny                         = std::stoi(m0.find(gdcm::Tag(0x0028, 0x0010))->second);
                nx                         = std::stoi(m0.find(gdcm::Tag(0x0028, 0x0011))->second);
                rect3d<int16_t, R>::dim_.x = nx;
                rect3d<int16_t, R>::dim_.y = ny;
                std::string  pixel_spacing(m0.find(gdcm::Tag(0x0028, 0x0030))->second);
                unsigned int deli = pixel_spacing.find_first_of('\\');
                dx_               = std::stod(pixel_spacing.substr(0, deli));
                dy_               = std::stod(pixel_spacing.substr(deli + 1));
            }

            ///< A map to search file path upon instance UID
            uid2file_.insert(std::make_pair(m0.find(gdcm::Tag(0x0008, 0x0018))->second, files_[i]));
        }

        rect3d<int16_t, R>::x_ = new R[nx];
        for (size_t i = 0; i < nx; ++i) {
            rect3d<int16_t, R>::x_[i] = x0 + dx_ * i;
        }
        rect3d<int16_t, R>::y_ = new R[nx];
        for (size_t i = 0; i < ny; ++i) {
            rect3d<int16_t, R>::y_[i] = y0 + dy_ * i;
        }
        std::cout << "Reading DICOM directory.. : Getting patient CT info.." << std::endl;
        std::cout << "Reading DICOM directory.. : Detected CT info is as below." << std::endl;
        std::cout << "CT (nx, ny, nz) -> (" << rect3d<int16_t, R>::dim_.x << ", "
                  << rect3d<int16_t, R>::dim_.y << ", " << rect3d<int16_t, R>::dim_.z << ")"
                  << std::endl;
        std::cout << "CT (dx, dy, dz) -> (" << dx_ << ", " << dy_ << ", "
                  << rect3d<int16_t, R>::z_[1] - rect3d<int16_t, R>::z_[0] << ")" << std::endl;
        std::cout << "CT (x, y, z) -> (" << rect3d<int16_t, R>::x_[0] << ", "
                  << rect3d<int16_t, R>::y_[0] << ", " << rect3d<int16_t, R>::z_[0] << ")"
                  << std::endl;
    }

    /// \brief Loads the pixel data from the DICOM files into the grid.
    ///
    /// This method reads the pixel data for each slice and populates the `data_`
    /// member of the base `rect3d` class. It also applies the rescale slope
    /// and intercept found in the DICOM metadata.
    CUDA_HOST
    virtual void
    load_data() {
        std::cout << "Reading DICOM directory.. : Loading patient CT pixel data.." << std::endl;
        size_t nb_voxels_2d = rect3d<int16_t, R>::dim_.x * rect3d<int16_t, R>::dim_.y;
        size_t nb_voxels_3d = nb_voxels_2d * rect3d<int16_t, R>::dim_.z;

        rect3d<int16_t, R>::data_.resize(nb_voxels_3d);
        float intercept = 0;
        float slope     = 1;

        for (size_t i = 0; i < rect3d<int16_t, R>::dim_.z; ++i) {

            gdcm::ImageReader reader;
            reader.SetFileName(files_[i].c_str());
            reader.Read();
            const gdcm::Image& img = reader.GetImage();

            //n_x * n_y * bytes = img.GetBufferLength()
            intercept = float(img.GetIntercept());
            slope     = float(img.GetSlope());

            gdcm::PixelFormat pixeltype = img.GetPixelFormat();

            switch (pixeltype) {
            case gdcm::PixelFormat::INT16: {
                img.GetBuffer((char*) &rect3d<int16_t, R>::data_[i * nb_voxels_2d]);
            } break;
            default:
                assert(0);
            }   //switch

        }   //for
        rect3d<int16_t, R>::data_ = rect3d<int16_t, R>::data_ * int16_t(slope) + intercept;
        std::cout << "Reading DICOM directory.. : Patient CT pixel data successfully loaded." << std::endl;
    }

    /// \brief Finds the x-index corresponding to a given x-coordinate.
    /// \param[in] x The x-coordinate.
    /// \return The index of the voxel in the x-dimension.
    inline virtual size_t
    find_c000_x_index(const R& x) {
        assert((x >= rect3d<int16_t, R>::x_[0]) &&
               (x < rect3d<int16_t, R>::x_[rect3d<int16_t, R>::dim_.x - 1]));
        assert(dx_ > 0);
        return floor(x / dx_);
    }

    /// \brief Finds the y-index corresponding to a given y-coordinate.
    /// \param[in] y The y-coordinate.
    /// \return The index of the voxel in the y-dimension.
    inline virtual size_t
    find_c000_y_index(const R& y) {
        assert((y >= rect3d<int16_t, R>::y_[0]) &&
               (y < rect3d<int16_t, R>::y_[rect3d<int16_t, R>::dim_.y - 1]));
        assert(dy_ > 0);
        return floor(y / dy_);
    }

    /// \brief Gets the pixel spacing in the x-dimension.
    /// \return The pixel spacing in x.
    inline virtual R
    get_dx() {
        return dx_;
    }

    /// \brief Gets the pixel spacing in the y-dimension.
    /// \return The pixel spacing in y.
    inline virtual R
    get_dy() {
        return dy_;
    }

    /// \brief Gets the array of pixel spacings in the z-dimension.
    /// \return A pointer to the array of z-spacings.
    inline virtual R*
    get_dz() {
        return dz_;
    }

    /// \brief Sets the pixel spacing in the x-dimension.
    /// \param[in] dx The new pixel spacing in x.
    inline virtual void
    set_dx(R dx) {
        dx_ = dx;
    }

    /// \brief Sets the pixel spacing in the y-dimension.
    /// \param[in] dy The new pixel spacing in y.
    inline virtual void
    set_dy(R dy) {
        dy_ = dy;
    }

    /// \brief Sets the array of pixel spacings in the z-dimension.
    /// \param[in] dz A pointer to the new array of z-spacings.
    inline virtual void
    set_dz(R* dz) {
        dz_ = dz;
    }

    /// \brief Sets the data (incorrectly named, should likely not be used).
    /// \param[in] dz A pointer to a data array.
    inline virtual void
    set_data(R* dz) {
        dz_ = dz;
    }

    /// \brief A friend function to clone the grid structure of a ct object.
    /// \tparam R0 The coordinate type of the source ct grid.
    /// \tparam T1 The pixel type of the destination grid.
    /// \tparam R1 The coordinate type of the destination grid.
    /// \param[in] src The source ct object whose structure is to be copied.
    /// \param[out] dest The destination `rect3d` object that will receive the structure.
    template<typename R0, typename T1, typename R1>
    friend void
    clone_ct_structure(ct<R0>& src, rect3d<T1, R1>& dest);
};

}   // namespace mqi

#endif
