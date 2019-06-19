//----------------------------------*-C++-*--------------------------------//
/*!
 * \file   nuclear_data/Nuclear_Data.cc
 * \author B. R. Ryan
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for nuclear_data/Nuclear_Data class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Nuclear_Data.hh"

namespace rtt_nuclear_data {

//----------------------------------------------------------------------------//
/*!
 * \brief Parses the cell_data block data from the mesh file via calls to
 *        private member functions.
 * \param meshfile Mesh file name.
 */
Nuclear_Data::Nuclear_Data(string const &filename) {
  ifstream ACEfile(filename.c_str(), std::ios::in);
  if (!ACEfile) {
    std::ostringstream buffer;
    buffer << "File " << filename << " could not be opened\n";
    throw std::invalid_argument(buffer.str());
  }
}
} // end namespace rtt_nuclear_data
