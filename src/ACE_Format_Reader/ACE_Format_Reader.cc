//----------------------------------*-C++-*--------------------------------//
/*!
 * \file   ACE_Format_Reader/ACE_Format_Reader.cc
 * \author B. R. Ryan
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for ACE_Format_Reader library.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ACE_Format_Reader.hh"
#include "ds++/Assert.hh"
#include <sstream>

namespace rtt_ACE_Format_reader {
/*!
 * \brief Constructs an ACE_Format_Reader object and parses ACE file
 * \param ACE_File Nuclear data file name.
 */
ACE_Format_Reader::ACE_Format_Reader(string const &ACE_File) {
  readFile(ACE_File);
};
/*!
 * \brief Parses the nuclear data file, stores in class data members
 * \param ACE_File Nuclear data file name.
 * \throw invalid_argument If the file cannot be opened
 */
void ACE_Format_Reader::readFile(const string &ACE_File) {
  const char *filename = ACE_File.c_str();
  ifstream acefile(filename, std::ios::in);
  if (!acefile) {
    std::ostringstream buffer;
    buffer << "File " << ACE_File << " could not be opened\n";
    throw std::invalid_argument(buffer.str());
  }

  try {
  } catch (rtt_dsxx::assertion &as) {
    std::cout << "Assertion thrown: " << as.what() << std::endl;
    Insist(false, as.what());
  }
}

} // end namespace rtt_ACE_Format_Reader

//---------------------------------------------------------------------------//
// end of ACE_Format_Reader/ACE_Format_Reader.cc
//---------------------------------------------------------------------------//
