//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/IpcressOdfmgOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 22 15:24210 2001
 * \brief  IpcressOdfmgOpacity templated class implementation file.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "IpcressOdfmgOpacity.hh"
#include "IpcressDataTable.hh"
#include "IpcressFile.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace rtt_cdi_ipcress {

// ------------ //
// Constructors //
// ------------ //

/*!
 * \brief Constructor for IpcressOdfmgOpacity object.
 *
 * See IpcressOdfmgOpacity.hh for details.
 */
IpcressOdfmgOpacity::IpcressOdfmgOpacity(
    std::shared_ptr<const IpcressFile> const &in_spIpcressFile,
    size_t in_materialID, rtt_cdi::Model in_opacityModel,
    rtt_cdi::Reaction in_opacityReaction, size_t numBands)
    : spIpcressFile(in_spIpcressFile), materialID(in_materialID),
      //numKeys( 0 ),
      vKnownKeys(), opacityModel(in_opacityModel),
      opacityReaction(in_opacityReaction), energyPolicyDescriptor("odfmg"),
      spIpcressDataTable(), groupBoundaries(), bandBoundaries(),
      reverseBands(false) {
  // Verify that the requested material ID is available in the specified
  // IPCRESS file.
  Insist(spIpcressFile->materialFound(materialID),
         std::string("The requested material ID is not available in the ") +
             std::string("specified Ipcress file."));

  // Retrieve keys available fo this material from the IPCRESS
  // file.
  vKnownKeys = spIpcressFile->listDataFieldNames(materialID);
  Check(vKnownKeys.size() > 0);

  // Create the data table object and fill it with the table data from the
  // IPCRESS file.  As far as the data table is concerned, this is MG data.
  spIpcressDataTable.reset(new IpcressDataTable("mg", opacityModel,
                                                opacityReaction, vKnownKeys,
                                                materialID, spIpcressFile));

  //load the group and band structure
  loadGroupsAndBands(numBands);

} // end of IpcressData constructor

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor for IpcressOdfmgOpacity object.
 *
 * See IpcressOdfmgOpacity.hh for details.
 */
IpcressOdfmgOpacity::IpcressOdfmgOpacity(std::vector<char> const &packed)
    : spIpcressFile(), materialID(0),
      //numKeys(0),
      vKnownKeys(), opacityModel(), opacityReaction(),
      energyPolicyDescriptor("odfmg"), spIpcressDataTable(), groupBoundaries(),
      bandBoundaries(), reverseBands(false) {
  Require(packed.size() >= 6 * sizeof(int));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack and check the descriptor
  int packed_descriptor_size = 0;
  unpacker >> packed_descriptor_size;
  Check(packed_descriptor_size > 0);

  // make a vector<char> for the packed descriptor
  std::vector<char> packed_descriptor(packed_descriptor_size);

  // unpack it
  std::string descriptor;
  for (size_t i = 0; i < static_cast<size_t>(packed_descriptor_size); i++)
    unpacker >> packed_descriptor[i];
  rtt_dsxx::unpack_data(descriptor, packed_descriptor);

  // make sure it is "odfmg"
  Insist(descriptor == "odfmg",
         "Tried to unpack a non-odfmg opacity in IpcressOdfmgOpacity.");

  // unpack the size of the packed filename
  int packed_filename_size = 0;
  unpacker >> packed_filename_size;

  // make a vector<char> for the packed filename
  std::vector<char> packed_filename(packed_filename_size);

  // unpack it
  std::string filename;
  for (size_t i = 0; i < static_cast<size_t>(packed_filename_size); i++)
    unpacker >> packed_filename[i];
  rtt_dsxx::unpack_data(filename, packed_filename);

  // unpack the material id
  int matid(0);
  unpacker >> matid;
  materialID = static_cast<size_t>(matid);

  // unpack the model and reaction
  int model = 0;
  int reaction = 0;
  unpacker >> model >> reaction;

  opacityModel = static_cast<rtt_cdi::Model>(model);
  opacityReaction = static_cast<rtt_cdi::Reaction>(reaction);

  // Load the number of bands
  int numBands = 0;
  unpacker >> numBands;
  Check(numBands > 0);

  Ensure(unpacker.get_ptr() == &packed[0] + packed.size());

  // build a new IpcressFile
  spIpcressFile.reset(new IpcressFile(filename));
  Check(spIpcressFile);

  // Verify that the requested material ID is available in the specified
  // IPCRESS file.
  Insist(spIpcressFile->materialFound(materialID),
         "Requested material ID is not found in the specified Ipcress file.");

  // Retrieve keys available fo this material from the IPCRESS file.
  vKnownKeys = spIpcressFile->listDataFieldNames(materialID);
  Check(vKnownKeys.size() > 0);

  // Create the data table object and fill it with the table data from the
  // IPCRESS file.  To the data table, this is a 'mg' object.
  spIpcressDataTable.reset(new IpcressDataTable("mg", opacityModel,
                                                opacityReaction, vKnownKeys,
                                                materialID, spIpcressFile));

  Ensure(spIpcressFile);
  Ensure(spIpcressDataTable);

  // set up the groups and bands
  loadGroupsAndBands(numBands);
}

//---------------------------------------------------------------------------//
void IpcressOdfmgOpacity::loadGroupsAndBands(size_t const numBands) {
  using rtt_dsxx::soft_equiv;
  using std::string;

  Require(numBands > 0);

  // we'll need to access the original boundaries a lot
  std::vector<double> const &origBoundaries =
      spIpcressDataTable->getGroupBoundaries();

  // number of groups*bands stored in the IPCRESS file
  size_t const numEffGroups = origBoundaries.size() - 1;

  // integer division here!
  size_t const numGroups = numEffGroups / numBands;

  // make sure that the groups and bands divide evenly
  Insist((numGroups * numBands) == (origBoundaries.size() - 1),
         string("The number of bands does not divide evenly into the ") +
             string("number of 'groups' stored in the IPCRESS file!"));

  // first, load the actual group structure
  groupBoundaries.assign(numGroups + 1, 0);
  for (size_t group = 0; group < numGroups + 1; group++) {
    Check(group * numBands < origBoundaries.size());
    groupBoundaries[group] = origBoundaries[group * numBands];
  }

  // next, check the band structure inside every group,
  // as given by tops (from highest opacity to lowest)
  for (size_t group = 0; group < numGroups; group++) {
    std::vector<double> currentBandStructure(numBands + 1);

    double const groupWidth =
        groupBoundaries[group + 1] - groupBoundaries[group];

    // figure out the band structure inside each group
    for (size_t band = 0; band <= numBands; band++) {
      Check(group * numBands + band < origBoundaries.size());
      // bands are stored between groups, really
      currentBandStructure[band] = origBoundaries[group * numBands + band];

      // subtract beginning group value, divide by width
      // (undo the intra-group interpolation that Tops does)
      currentBandStructure[band] =
          (currentBandStructure[band] - groupBoundaries[group]) / groupWidth;
    }

    if (group == 0) {
      // the first group will give us the band structure
      bandBoundaries = currentBandStructure;
    } else {
      // the rest of the bands should darned well be identical; check them.
      // Also, loosen the tolerance for more groups.
      if (!soft_equiv(bandBoundaries.begin(), bandBoundaries.end(),
                      currentBandStructure.begin(), currentBandStructure.end(),
                      numEffGroups * 1.e-12)) {
        std::cerr << "Band boundaries do not match.\n"
                  << "First (reverse, tops-style) band structure, " << group
                  << "th band structure: \n";
        for (size_t i = 0; i <= numBands; i++) {
          std::cerr << i + 1 << "\t" << std::setprecision(14)
                    << bandBoundaries[i] << "\t" << std::setprecision(14)
                    << currentBandStructure[i] << std::endl;
        }

        throw rtt_dsxx::assertion(
            string("Band boundaries do not match.  ") +
            string("IpcressOdfmgOpacity::loadGroupsAndBands"));
        ;
      }
    }
  } // end internal band checking and assignment

  Check(bandBoundaries.size() == numBands + 1);

  // check the opacities at some temperature and density on the grid
  // to find whether the order is increasing or decreasing in each group
  if (numBands > 1) {
    // pick some points in the middle of the data, probably not necessary
    size_t groupPoint = numGroups / 2;
    size_t temperaturePoint = getNumTemperatures() / 2;
    size_t densityPoint = getNumDensities() / 2;

    std::vector<std::vector<double>> opacity;
    opacity = getOpacity(getTemperatureGrid()[temperaturePoint],
                         getDensityGrid()[densityPoint]);

    Check(reverseBands == false);

    if (opacity[groupPoint][numBands - 1] < opacity[groupPoint][0])
      reverseBands = true;
  } else // don't bother reversing if only one band
  {
    reverseBands = false;
  }

  if (reverseBands) {
    // now reverse the band structure from tops-style to traditional ODF
    // (i.e., highest opacity occupies the highest band)
    std::vector<double> bandWidths(numBands);

    // original band widths
    for (size_t band = 0; band < numBands; band++) {
      bandWidths[band] = bandBoundaries[band + 1] - bandBoundaries[band];
    }
    // reverse their order
    std::reverse(bandWidths.begin(), bandWidths.end());

    // re-accumulate into boundaries
    bandBoundaries[0] = 0.0;
    for (size_t band = 0; band < numBands; band++) {
      bandBoundaries[band + 1] = bandBoundaries[band] + bandWidths[band];
    }
  }

  if (soft_equiv(bandBoundaries[0], 0.0))
    bandBoundaries[0] = 0.0;
  else
    throw rtt_dsxx::assertion(
        string("IpcressOdfmgOpacity::loadGroupsAndBands :: ") +
        string("Bad bad bad band boundaries! ") +
        string("Lowest band boundary is non-zero."));

  if (soft_equiv(bandBoundaries[numBands], 1.0))
    bandBoundaries[numBands] = 1.0;
  else
    throw rtt_dsxx::assertion(
        string("IpcressOdfmgOpacity::loadGroupsAndBands :: ") +
        string("Bad bad bad band boundaries! Highest band boundary ") +
        string("should be one."));
}

/*!
 * \brief Default IpcressOpacity() destructor.
 *
 * \sa This is required to correctly release memory when a IpcressOdfmgOpacity
 *     is destroyed.  This constructor's definition must be declared in the
 *     implementation file so that we can avoid including too many header
 *     files
 */
IpcressOdfmgOpacity::~IpcressOdfmgOpacity() {
  // empty
}

// --------- //
// Accessors //
// --------- //

/*!
 * \brief Returns a "plain English" description of the opacity data that this
 *     class references. (e.g. "Odfmg Rosseland Scattering".)
 *
 * The definition of this function is not included here to prevent the
 *     inclusion of the IpcressFile.hh definitions within this header file.
 */
std::string IpcressOdfmgOpacity::getDataDescriptor() const {
  // call the correct function from the IpcressDataTable object.
  return spIpcressDataTable->getDataDescriptor();
}

/*!
 * \brief Returns the name of the associated IPCRESS file.
 *
 * The definition of this function is not included here to prevent the
 *     inclusion of the IpcressFile.hh definitions within this header file.
 */
std::string IpcressOdfmgOpacity::getDataFilename() const {
  return spIpcressFile->getDataFilename();
}

/*!
 * \brief Opacity accessor that returns a 2-D vector (bands * groups) of
 *     opacities that corresponds to the provided temperature and density.
 */
std::vector<std::vector<double>>
IpcressOdfmgOpacity::getOpacity(double targetTemperature,
                                double targetDensity) const {
  // number of groups*bands in this multigroup set (effective groups).
  size_t const numEffGroups = spIpcressDataTable->getNumGroupBoundaries() - 1;
  size_t const numBands = getNumBands();
  size_t const numGroups = getNumGroups();

  Check(numEffGroups == numGroups * numBands);

  // temporary opacity vector used by the wrapper: TOPS stores the
  // group*band information in a flat list, with the band opacities in a
  // group in reverse order
  std::vector<double> tempOpacity(numEffGroups);

  // logarithmic interpolation:
  for (size_t g = 0; g < numEffGroups; ++g) {
    tempOpacity[g] =
        spIpcressDataTable->interpOpac(targetTemperature, targetDensity, g);
    Check(tempOpacity[g] >= 0.0);
  }
  // logarithmic interpolation:
  // tempOpacity = wrapper::wgintmglog(
  //     spIpcressDataTable->getLogTemperatures(),
  //     spIpcressDataTable->getNumTemperatures(),
  //     spIpcressDataTable->getLogDensities(),
  //     spIpcressDataTable->getNumDensities(),
  //     spIpcressDataTable->getNumGroupBoundaries(),
  //     spIpcressDataTable->getLogOpacities(),
  //     spIpcressDataTable->getNumOpacities(),
  //     std::log( targetTemperature ),
  //     std::log( targetDensity ) );

  // groups * bands opacity 2D-vector to be returned by this function
  std::vector<std::vector<double>> opacity(numGroups);

  for (size_t group = 0; group < numGroups; group++) {
    // set the size of each group array to the number of bands
    opacity[group].resize(numBands);

    // fill each band opacity, in reverse order
    for (size_t band = 0; band < numBands; band++) {
      if (reverseBands) {
        opacity[group][band] =
            tempOpacity[group * numBands + (numBands - band - 1)];
      } else {
        opacity[group][band] = tempOpacity[group * numBands + band];
      }
    }
  }

  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband opacity 2-D
 *     vectors that correspond to the provided vector of temperatures and a
 *     single density value.
 */
std::vector<std::vector<std::vector<double>>>
IpcressOdfmgOpacity::getOpacity(std::vector<double> const &targetTemperature,
                                double targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(
      targetTemperature.size());

  for (size_t i = 0; i < targetTemperature.size(); ++i)
    opacity[i] = getOpacity(targetTemperature[i], targetDensity);
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided
 *     temperature and a vector of density values.
 */
std::vector<std::vector<std::vector<double>>> IpcressOdfmgOpacity::getOpacity(
    double targetTemperature, const std::vector<double> &targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(targetDensity.size());

  //call our regular getOpacity function for every target density
  for (size_t i = 0; i < targetDensity.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a vector of temperatures that define the cached opacity data
 *     table.
 */
std::vector<double> IpcressOdfmgOpacity::getTemperatureGrid() const {
  return spIpcressDataTable->getTemperatures();
}

//---------------------------------------------------------------------------//
//! Returns the size of the temperature grid.
size_t IpcressOdfmgOpacity::getNumTemperatures() const {
  return spIpcressDataTable->getNumTemperatures();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns a vector of densities that define the cached opacity data
 *     table.
 */
std::vector<double> IpcressOdfmgOpacity::getDensityGrid() const {
  return spIpcressDataTable->getDensities();
}

//---------------------------------------------------------------------------//
//!  Returns the size of the density grid.
size_t IpcressOdfmgOpacity::getNumDensities() const {
  return spIpcressDataTable->getNumDensities();
}

// ------- //
// Packing //
// ------- //

//---------------------------------------------------------------------------//
/*!
 * Pack the IpcressOdfmgOpacity state into a char string represented by a
 * vector<char>. This can be used for persistence, communication, etc. by
 * accessing the char * under the vector (required by implication by the
 * standard) with the syntax &char_string[0]. Note, it is unsafe to use
 * iterators because they are \b not required to be char *.
 */
std::vector<char> IpcressOdfmgOpacity::pack() const {
  using std::vector;
  using std::string;

  Require(spIpcressFile);

  // pack up the energy policy descriptor
  vector<char> packed_descriptor;
  rtt_dsxx::pack_data(energyPolicyDescriptor, packed_descriptor);

  // pack up the ipcress file name
  string filename = spIpcressFile->getDataFilename();
  vector<char> packed_filename;
  rtt_dsxx::pack_data(filename, packed_filename);

  // determine the total size: 3 ints (reaction, model, material id) + 2
  // ints for packed_filename size and packed_descriptor size + 1 int for
  // number of bands + char in packed_filename and packed_descriptor
  size_t size =
      6 * sizeof(int) + packed_filename.size() + packed_descriptor.size();

  // make a container to hold packed data
  vector<char> packed(size);

  // make a packer and set it
  rtt_dsxx::Packer packer;
  packer.set_buffer(size, &packed[0]);

  // pack the descriptor
  packer << static_cast<int>(packed_descriptor.size());
  for (size_t i = 0; i < packed_descriptor.size(); i++)
    packer << packed_descriptor[i];

  // pack the filename (size and elements)
  packer << static_cast<int>(packed_filename.size());
  for (size_t i = 0; i < packed_filename.size(); i++)
    packer << packed_filename[i];

  // pack the material id
  packer << static_cast<int>(materialID);

  // pack the model and reaction
  packer << static_cast<int>(opacityModel) << static_cast<int>(opacityReaction);

  //pack the number of bands
  packer << static_cast<int>(getNumBands());

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

} // end namespace rtt_cdi_ipcress

//---------------------------------------------------------------------------//
// end of IpcressOdfmgOpacity.cc
//---------------------------------------------------------------------------//
