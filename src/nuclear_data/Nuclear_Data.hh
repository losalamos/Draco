//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nuclear_data/Nuclear_Data.hh
 * \author B. R. Ryan
 * \date   Wed Jun 19 2019
 * \brief  Header file for Nuclear_Data class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __nuclear_data_Nuclear_Data_hh__
#define __nuclear_data_Nuclear_Data_hh__

namespace rtt_nuclear_data {
//===========================================================================//
/*!
 * \class Nuclear_Data
 *
 * \brief A class for storing the data from one ACE nuclear data file.
 *
 *\sa The Nuclear_Data class contains data members for storing all the data in
 *    an ACE nuclear data file. The name of an ACE file is passed to the class
 *    constructor, which then reads that file and populates the class's data
 *    members. Public functions are provided for accessing this data.
 */
//===========================================================================//
class DLL_PUBLIC_Nuclear_Data Nuclear_Data {
public:
  Nuclear_Data(string const &filename);
  ~NuclearData() {}
};
} // end namespace rtt_nuclear_data
