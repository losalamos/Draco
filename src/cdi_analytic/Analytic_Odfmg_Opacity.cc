//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Analytic_Odfmg_Opacity class member definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Analytic_Odfmg_Opacity.hh"
#include "ds++/Packing_Utils.hh"

namespace rtt_cdi_analytic {

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic multigroup opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * The group structure (in keV) must be provided by the groups argument.  The
 * number of Analytic_Opacity_Model objects given in the models argument must
 * be equal to the number of groups.
 *
 * \param groups vector containing the group boundaries in keV from lowest to
 * highest
 *
 * \param models vector containing shared_ptrs to Analytic_Model derived types
 * for each group, the size should be groups.size() - 1
 *
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 *
 */
Analytic_Odfmg_Opacity::Analytic_Odfmg_Opacity(const sf_double &groups,
                                               const sf_double &bands,
                                               rtt_cdi::Reaction reaction_in,
                                               rtt_cdi::Model model_in)
    : groupBoundaries(groups), reaction(reaction_in), model(model_in),
      bandBoundaries(bands) {
  Require(reaction == rtt_cdi::TOTAL || reaction == rtt_cdi::ABSORPTION ||
          reaction == rtt_cdi::SCATTERING);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and Analytic_Odfmg_Opacity from a
 * vector<char> that was created by a call to pack().  It can only rebuild
 * Analytic_Model types that have been registered in the
 * rtt_cdi_analytic::Opacity_Models enumeration.
 */
Analytic_Odfmg_Opacity::Analytic_Odfmg_Opacity(const sf_char &packed)
    : groupBoundaries(0), reaction(), model(),
      bandBoundaries(std::vector<double>()) {
  // the packed size must be at least 5 integers (number of groups, number of
  // bands, reaction type, model type, analytic model indicator)
  Require(packed.size() >= 5 * sizeof(int));

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // register the unpacker
  unpacker.set_buffer(packed.size(), &packed[0]);

  // unpack the number of group boundaries
  int ngrp_bounds = 0;
  unpacker >> ngrp_bounds;
  // int num_groups  = ngrp_bounds - 1;

  // make the group boundaries and model vectors
  groupBoundaries.resize(ngrp_bounds);

  // unpack the group boundaries
  for (int i = 0; i < ngrp_bounds; i++)
    unpacker >> groupBoundaries[i];

  // unpack the number of band boundaries
  int nband_bounds = 0;
  unpacker >> nband_bounds;

  // make the group boundaries and model vectors
  bandBoundaries.resize(nband_bounds);

  // unpack the group boundaries
  for (int i = 0; i < nband_bounds; i++)
    unpacker >> bandBoundaries[i];

  // unpack the reaction and model type
  int reaction_int, model_int;
  unpacker >> reaction_int >> model_int;

  // assign the reaction and model type
  reaction = static_cast<rtt_cdi::Reaction>(reaction_int);
  model = static_cast<rtt_cdi::Model>(model_int);
  Check(reaction == rtt_cdi::ABSORPTION || reaction == rtt_cdi::SCATTERING ||
        reaction == rtt_cdi::TOTAL);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic odfmg opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char
 * array (represented by a vector<char>).  The Analytic_Opacity_Model derived
 * class must have a pack function; this is enforced by the virtual
 * Analytic_Opacity_Model base class.
 */
Analytic_Odfmg_Opacity::sf_char Analytic_Odfmg_Opacity::pack() const {
  // make a packer
  rtt_dsxx::Packer packer;

  // now add up the total size; number of groups + 1 size_t for number of
  // groups, number of bands + 1 size_t for number of
  // bands, number of models + size in each model + models, 1 size_t for
  // reaction type, 1 size_t for model type
  size_t size = 4 * sizeof(int) + groupBoundaries.size() * sizeof(double) +
                bandBoundaries.size() * sizeof(double);

  // make a char array
  sf_char packed(size);

  // set the buffer
  packer.set_buffer(size, &packed[0]);

  // pack the number of groups and group boundaries
  packer << static_cast<int>(groupBoundaries.size());
  for (size_t i = 0; i < groupBoundaries.size(); i++)
    packer << groupBoundaries[i];

  // pack the number of bands and band boundaries
  packer << static_cast<int>(bandBoundaries.size());
  for (size_t i = 0; i < bandBoundaries.size(); i++)
    packer << bandBoundaries[i];

  // now pack the reaction and model type
  packer << static_cast<int>(reaction) << static_cast<int>(model);

  Ensure(packer.get_ptr() == &packed[0] + size);
  return packed;
}

//---------------------------------------------------------------------------//
unsigned Analytic_Odfmg_Opacity::packed_size() const {
  // This must match the size calculated in the previous function
  return 4 * sizeof(int) + groupBoundaries.size() * sizeof(double) +
         bandBoundaries.size() * sizeof(double);
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Analytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
