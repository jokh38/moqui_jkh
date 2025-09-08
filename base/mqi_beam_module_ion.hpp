#ifndef MQI_BEAM_MODULE_ION_H
#define MQI_BEAM_MODULE_ION_H

/// \file mqi_beam_module_ion.hpp
///
/// \brief Interprets DICOM-RT Ion beam modules for plans and treatment records.
///
/// This file provides the `beam_module_ion` class, which is responsible for parsing
/// and managing data from DICOM-RT Ion Plan (RTI) and Ion Beam Treatment Record (RTIBTR)
/// modules. It extracts information about scan spots, energies, and weights.
///
/// \see http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.8.25.html for RTI
/// \see http://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.8.26.html for RTIBTR

#include <moqui/base/mqi_beam_module.hpp>

namespace mqi
{

/// \enum SCAN_MODE
/// \brief Defines the scan mode for an ion beam, based on DICOM tag (300A,0308).
typedef enum
{
    NONE           = 0, ///< No scan mode specified.
    UNIFORM        = 1, ///< Uniform scanning.
    MODULATED      = 2, ///< Modulated scanning.
    MODULATED_SPEC = 3  ///< Modulated scanning with specific parameters.
} SCAN_MODE;

/// \struct logfile_t
/// \brief Represents data for a single field in a treatment log file.
/// \note Added by Chanil Jeon of SMC (2023-11-02).
struct logfile_t
{
    std::vector<float> posX;    ///< X-positions of spots.
    std::vector<float> posY;    ///< Y-positions of spots.
    std::vector<int>   muCount; ///< Monitor units or particle counts for each spot.
};

/// \struct logfiles_t
/// \brief Represents data for all energy layers in a treatment log file.
/// \note Added by Chanil Jeon of SMC (2023-11-02).
struct logfiles_t
{
    std::vector<std::vector<float>>   beamEnergyInfo; ///< Information about beam energy for each layer.
    std::vector<std::vector<logfile_t>> beamInfo;       ///< Detailed spot information for each layer.
};

/// \class beam_module_ion
/// \brief A class for handling RT-ION beams from DICOM plans and treatment records.
///
/// This class inherits from `beam_module` and specializes in parsing the Ion Control
/// Point Sequence from DICOM-RT Ion objects. It extracts all necessary parameters for
/// simulating a pencil beam scanning treatment, including spot positions, energies,
/// FWHM, and meterset weights.
class beam_module_ion : public beam_module
{
public:
    /// \struct spot
    /// \brief A user-defined type representing a single scan spot.
    typedef struct {
        float e;        ///< The spot energy in MeV.
        float x;        ///< The x-position in mm.
        float y;        ///< The y-position in mm.
        float fwhm_x;   ///< The Full-Width at Half-Maximum in the x-direction.
        float fwhm_y;   ///< The Full-Width at Half-Maximum in the y-direction.
        float meterset; ///< The meterset weight (e.g., MU or number of particles).
    } spot;

    /// \struct logspot
    /// \brief Represents a single spot from a treatment log file for MC simulation.
    /// \note Added in 2023 for SMC Log-file based MC simulation.
    typedef struct {
        float e;       ///< The spot energy in MeV.
        float x;       ///< The x-position in mm.
        float y;       ///< The y-position in mm.
        float muCount; ///< The particle count for Monte Carlo simulation.
    } logspot;

protected:
    /// \brief A vector containing the number of spots in each energy layer.
    std::vector<int> nb_spots_per_layer_;

    /// \brief A vector containing all spots in the order they are delivered.
    std::vector<spot> sequence_;

    /// \brief The name of the beam model or tune ID.
    std::string tune_id_;

public:
    /// \brief Constructs a beam module for RT-Ion.
    ///
    /// This constructor parses the provided DICOM dataset to extract the ion beam
    /// parameters. It handles both plan (IONPLAN) and treatment record (IONRECORD)
    /// modalities.
    ///
    /// \param d A pointer to the DICOM dataset, which should be an item from the
    ///          IonBeamSequence for a plan (C.8.8.25-1) or a record (C.8.8.26-1).
    /// \param m The modality type, which must be either `IONPLAN` or `IONRECORD`.
    beam_module_ion(const mqi::dataset* d, mqi::modality_type m) : beam_module(d, m) {
        /// Initializes containers
        std::vector<int>         nb_pts(1);
        std::vector<float>       energy(1);
        std::vector<float>       fwhm_xy(2);
        std::vector<std::string> tune_id(1);
        std::vector<float>       xy;
        std::vector<float>       weight;

        int layer_nb = 0;

        std::string str_weight;

        /// Checks modality type
        switch (modality_) {
        case mqi::modality_type::IONPLAN:
            str_weight = "ScanSpotMetersetWeights";
            break;
        case mqi::modality_type::IONRECORD:
            str_weight = "ScanSpotMetersetsDelivered";
            break;
        default:
            throw std::runtime_error("Wrong ION type");
        }

        /// Fill the containers from the DICOM dataset
        /// As each layer consists of pair of ion-controls and even-layers have all zero weights.
        /// So we drop even layer.
        auto seq_tags = &mqi::seqtags_per_modality.at(m);
        auto ictrl    = (*ds_)(seq_tags->at("ctrl"));
        for (auto b : ictrl) {

            if ((layer_nb++) % 2) continue;

            b->get_values("ScanSpotTuneID", tune_id);
            b->get_values("NominalBeamEnergy", energy);
            b->get_values("NumberOfScanSpotPositions", nb_pts);
            b->get_values("ScanningSpotSize", fwhm_xy);
            b->get_values("ScanSpotPositionMap", xy);
            b->get_values(str_weight.c_str(), weight);

            for (int j = 0; j < nb_pts[0]; ++j) {
                tune_id_ = tune_id[0];
                sequence_.push_back(
                  { energy[0], xy[j * 2], xy[j * 2 + 1], fwhm_xy[0], fwhm_xy[1], weight[j] });
            }   //per spot
            nb_spots_per_layer_.push_back(nb_pts[0]);
        }   //per layer
    }

    /// \brief Destroys the beam_module_ion object.
    ~beam_module_ion() {
        ;
    }

    /// \brief Returns a constant pointer to the vector of spot counts per layer.
    ///
    /// \return A pointer to a vector where each element is the number of spots
    ///         in the corresponding energy layer.
    const std::vector<int>*
    get_nb_spots_per_layer(void) const {
        return &nb_spots_per_layer_;
    }

    /// \brief Returns a constant pointer to the spot sequence vector.
    ///
    /// \return A pointer to a vector containing all the spots in delivery order.
    const std::vector<spot>*
    get_sequence(void) const {
        return &sequence_;
    }

    /// \brief Returns the tune ID of the beam model.
    ///
    /// \return A constant string containing the tune ID.
    const std::string
    get_tune_id(void) const {
        return tune_id_;
    }

    /// \brief Prints the details of the spot sequence to the console.
    ///
    /// This method iterates through the spot sequence and prints the energy, position,
    /// FWHM, and meterset weight for each spot. Useful for debugging.
    void
    dump() const {
        std::cout << "dump:spotmap, size:" << sequence_.size() << std::endl;
        for (auto i : sequence_) {
            std::cout << "spot (E,X,Y,Sx,Sy,W): " << i.e << ", " << i.x << ", " << i.y << ", "
                      << i.fwhm_x << ", " << i.fwhm_y << ", " << i.meterset << std::endl;
        }
    }
};

}   // namespace mqi

#endif
