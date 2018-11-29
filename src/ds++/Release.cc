//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Release.cc
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Provides the function definition for Release.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Release.hh"
#include "ds++/config.h"
#include <cstring> // memcpy
#include <sstream>

namespace rtt_dsxx {

//---------------------------------------------------------------------------//
/*!
 * \brief Format list of authors to do correct line breaks and ensures total
 *        line length is less than a specified maximum.
 *
 * \arg[in] maxlen Maximum line length
 * \arg[in] line_name Category title
 * \arg[in] list of developers
 * \return A formatted message.
 */
std::string print_devs(size_t const maxlinelen, std::string const &line_name,
                       mmdevs const &devs) {
  // indent subsequent lines by this many spaces.
  size_t const indent(5);
  std::string current_line(line_name);

  // temporary storage
  std::ostringstream msg;

  for (mmdevs::const_iterator it = devs.begin(); it != devs.end();) {
    std::string const name = it->second;
    if (current_line.length() + name.length() + 2 > maxlinelen) {
      // flush current line to the real output
      msg << current_line << std::endl;
      // reset the string that contains the current line.
      current_line.clear();
      // tab over to the colon
      current_line = std::string(indent, ' ');
    }
    // add the current developer to the list.
    if (++it == devs.end())
      current_line += std::string("and ") + name;
    else
      current_line += name + std::string(", ");
  }
  msg << current_line << "." << std::endl;

  return msg.str();
}

//---------------------------------------------------------------------------//
// function definition for Release, define the local version number for this
// library in the form ds_#.#.# in pkg_version variable
const std::string release() {
  std::ostringstream pkg_release;
  // Name and version
  pkg_release << "Draco-" << Draco_VERSION_MAJOR << "_" << Draco_VERSION_MINOR
              << "_" << Draco_VERSION_PATCH;

  // build date and type
  std::string const build_date(Draco_BUILD_DATE);
  std::string const build_type(BUILD_TYPE);
  pkg_release << ", build date " << build_date << "; build type: " << build_type
#ifdef DBC
              << "; DBC: " << DBC
#endif
#ifdef DRACO_DIAGNOSTICS
              << "; DRACO_DIAGNOSTICS: " << DRACO_DIAGNOSTICS
#endif
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#ifdef FPETRAP_SUPPORTED
              << "; FPE_TRAP: ON"
#endif
#endif
      ;

  return pkg_release.str();
}

//---------------------------------------------------------------------------//
/*! \brief Return a list of Draco contributing authors
 *
 * \param[in] use_doxygen_formatting If true, use extra decoration in the
 *              output.
 *
 * Data is collected from git (see regression/alist.sh) based on LOC
 * added/removed. Because the git repository only includes code provided
 * starting at draco-6_0_0, all LOC were attributed to KT at draco-6_0_0 since
 * he converted the svn repo to git. The remaining numbers are computed by
 * couting LOC added/removed since draco-6_0_0.
 */
const std::string author_list(bool const use_doxygen_formatting) {
  std::stringstream alist;

  mmdevs current_developers;
  // not totally fair... KT got credit for LOC when svn repository was converted
  // to git.
  current_developers.insert(fomdev(227428, "Kelly G. Thompson"));
  current_developers.insert(fomdev(12535, "Kent G. Budge"));
  current_developers.insert(fomdev(3464, "Ryan T. Wollaeger"));
  current_developers.insert(fomdev(3070, "James S. Warsa"));
  current_developers.insert(fomdev(2627, "Alex R. Long"));
  current_developers.insert(fomdev(1000, "Kendra P. Keady"));
  current_developers.insert(fomdev(398, "Jae H. Chang"));
  current_developers.insert(fomdev(244, "Matt A. Cleveland"));
  current_developers.insert(fomdev(120, "Andrew T. Till"));
  current_developers.insert(fomdev(108, "Tim Kelley"));

  mmdevs prior_developers;

  prior_developers.insert(fomdev(4868, "Jeff D. Densmore"));
  prior_developers.insert(fomdev(4058, "Gabriel M. Rockefeller"));
  prior_developers.insert(fomdev(2289, "Allan B. Wollaber"));
  prior_developers.insert(fomdev(1450, "Rob B. Lowrie"));
  prior_developers.insert(fomdev(995, "Lori A. Pritchett-Sheats"));
  prior_developers.insert(fomdev(313, "Paul W. Talbot"));
  prior_developers.insert(fomdev(262, "Katherine J. Wang"));
  prior_developers.insert(fomdev(78, "Peter Ahrens"));
  prior_developers.insert(fomdev(25, "Daniel Holladay"));
  prior_developers.insert(fomdev(9, "Massimiliano Rosa"));
  prior_developers.insert(fomdev(7, "Todd J. Urbatsch"));

  // Previous authors with no current LOC attribution:
  prior_developers.insert(fomdev(1, "Gabe Rockefeller"));
  prior_developers.insert(fomdev(1, "Jeff Furnish"));
  prior_developers.insert(fomdev(1, "John McGhee"));
  prior_developers.insert(fomdev(1, "Kris C. Garrett"));
  prior_developers.insert(fomdev(1, "Mike Buksas"));
  prior_developers.insert(fomdev(1, "Nick Myers"));
  prior_developers.insert(fomdev(1, "Paul Henning"));
  prior_developers.insert(fomdev(1, "Randy Roberts"));
  prior_developers.insert(fomdev(1, "Seth Johnson"));
  prior_developers.insert(fomdev(1, "Todd Adams"));
  prior_developers.insert(fomdev(1, "Tom Evans"));

  size_t maxlinelen(80);
  std::string line_name("CCS-2 Draco Team: ");

  if (use_doxygen_formatting) {
    maxlinelen = 400;
    alist << "\n\\par " << line_name << "\n\n";
    line_name = "";
  }

  alist << rtt_dsxx::print_devs(maxlinelen, line_name, current_developers);
  alist << "\n";

  line_name = std::string("Prior Contributers: ");
  if (use_doxygen_formatting) {
    alist << "\\par " << line_name << "\n\n";
    line_name = "";
  }
  alist << rtt_dsxx::print_devs(maxlinelen, line_name, prior_developers);

  return alist.str();
}

//---------------------------------------------------------------------------//
//! Print a Copyright note with an author list:
const std::string copyright() {
  std::ostringstream msg;

  msg << author_list() << "\n"
      << "Copyright (C) 2016-2018 Los Alamos National Security, LLC. "
         "(LA-CC-16-016)"
      << std::endl;

  return msg.str();
}

} // namespace rtt_dsxx

//---------------------------------------------------------------------------//
//! This version can be called by Fortran and wraps the C++ version.
extern "C" void ec_release(char *release_string, size_t maxlen) {
  std::string tmp_rel = rtt_dsxx::release();
  if (tmp_rel.size() >= static_cast<size_t>(maxlen)) {
    tmp_rel = tmp_rel.substr(0, maxlen - 1);
  }
  std::memcpy(release_string, tmp_rel.c_str(), tmp_rel.size() + 1);
  return;
}

//---------------------------------------------------------------------------//
// end of Release.cc
//---------------------------------------------------------------------------//
