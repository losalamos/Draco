//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/DracoTerminal.cc
 * \author Kelly Thompson
 * \date   Sat Oct 05 2019
 * \brief  Wrapper for a Terminal class that provides colored output.
 * \note   https://github.com/certik/terminal/blob/master/terminal.h
 *
 * \todo Consider an enum class for colors that derives from $LS_COLORS on
 *       Linux.  This would allow color selection based on users's terminal
 *       colors (e.g.: light vs dark scheme). */
//---------------------------------------------------------------------------//

#include "DracoTerminal.hh"

//----------------------------------------------------------------------------//
// Initialize pointer to zero so that it can be initialized in first call to
// getInstance
Term::DracoTerminal *Term::DracoTerminal::instance = nullptr;

//----------------------------------------------------------------------------//
//! Global bool to enable/disable color
bool Term::DracoTerminal::useColor = true;

//----------------------------------------------------------------------------//
/*! \brief Initialize other static const data (color map)
 *
 * \todo Eval ENV{GCC_COLORS} to override the defaults listed below? This
 *       probably would need to be done by cmake and saved to ds++/config.h.
 *
 * This should provide a mechanism to distinguish between light and dark
 * terminal backgrounds.
 *
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Message-Formatting-Options.html
 */
const std::vector<uint32_t> Term::DracoTerminal::error = {01, 31};
const std::vector<uint32_t> Term::DracoTerminal::warning = {01, 35};
const std::vector<uint32_t> Term::DracoTerminal::note = {01, 36};
const std::vector<uint32_t> Term::DracoTerminal::quote = {01};
const std::vector<uint32_t> Term::DracoTerminal::pass = {32};
const std::vector<uint32_t> Term::DracoTerminal::fail = {01, 31};
const std::vector<uint32_t> Term::DracoTerminal::reset = {00, 39};

//----------------------------------------------------------------------------//
//! specialization for std::vector<std::string>
std::string Term::ccolor(std::vector<uint32_t> const &value) {
  std::string retVal;
  if (Term::DracoTerminal::is_initialized() && Term::DracoTerminal::useColor)
    for (uint32_t const &it : value)
      retVal += "\033[" + std::to_string(static_cast<int>(it)) + "m";
  return retVal;
}

//---------------------------------------------------------------------------//
// end of ds++/DracoTerminal.cc
//---------------------------------------------------------------------------//
