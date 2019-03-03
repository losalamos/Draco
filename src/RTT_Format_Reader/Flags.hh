//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/Flags.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Flags class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_Flags_hh__
#define __RTT_Format_Reader_Flags_hh__

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {

//============================================================================//
/*!
 * \class Flags
 * \brief Controls parsing, storing, and accessing the data contained in the
 *        node, side, and cell flag blocks of the mesh file.
 */
//============================================================================//
class Flags {

  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<string> vector_str;
  typedef std::vector<int> vector_int;

  size_t nflags;
  string name;
  vector_int flag_nums;
  vector_str flag_names;

public:
  Flags(size_t nflags_, const string &name_)
      : nflags(nflags_), name(name_), flag_nums(nflags), flag_names(nflags) {}
  ~Flags() {}

  void readFlags(ifstream &meshfile);

public:
  /*!
   * \brief Validates the specified flag index.
   * \param flag Flag index.
   * \return The existance of the flag.
   */
  bool allowed_flag(size_t flag) const {
    return flag_nums.end() !=
           std::find(flag_nums.begin(), flag_nums.end(), flag);
  }

  /*!
   * \brief Returns the name of specified flag type
   * \return The Flag type name.
   */
  string getFlagType() const { return name; }

  /*!
   * \brief Returns the flag number associated with the specified flag index.
   * \param flag Flag index.
   * \return The Flag number.
   */
  int getFlagNumber(size_t flag) const { return flag_nums[flag]; }

  /*!
   * \brief Returns the flag name associated with the specified index.
   * \param flag Flag index.
   * \return The flag name.
   */
  string getFlagName(size_t flag) const { return flag_names[flag]; }

  /*!
   * \brief Returns the number of flags.
   * \return The number of flags.
   */
  size_t getFlagSize() const { return nflags; }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_Flags_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/Flags.hh
//---------------------------------------------------------------------------//
