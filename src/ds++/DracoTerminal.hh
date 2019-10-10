//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/DracoTerminal.hh
 * \author Kelly Thompson
 * \date   Sat Oct 05 2019
 * \brief  Wrapper for a Terminal class that provides colored output.
 * \note   https://github.com/certik/terminal/blob/master/terminal.h
 *
 * \todo Consider an enum class for colors that derives from $LS_COLORS on
 *       Linux.  This would allow color selection based on users's terminal
 *       colors (e.g.: light vs dark scheme). */
//---------------------------------------------------------------------------//

#ifndef dsxx_dracoterminal_hh
#define dsxx_dracoterminal_hh

#include "terminal/terminal.h"
#include "ds++/config.h"

namespace Term {

//===========================================================================//
/*!
 * \class DracoTerminal
 * \brief Global scope singleton object to ensure terminal setup/teardown is 
 *        done correctly.
 *
 * This will self construct on first call to Term::ccolor and it will remain
 * in scope until the program exits.
 */
class DracoTerminal {
  //! Private pointer to this object (this defines the singleton)
  DLL_PUBLIC_dsxx static DracoTerminal *instance;

  //! Private constructor so that no objects can be created.
  DracoTerminal() { /*empty*/
  }

  /*! \brief Construct a terminal object on creation.  This is will do some 
   *         terminal initialization that is required for some older software 
   *         (like Windows cmd.exe). */
  Term::Terminal term;

public:
  /*! \brief Calling this public function should always return true.
   *
   * If this singleton has not be created, then it will be created. Otherwise,
   * just check that the pointer to instance is valid and return true.
   */
  DLL_PUBLIC_dsxx static bool is_initialized() {
    if (instance == nullptr)
      instance = new DracoTerminal;
    return instance != nullptr;
  }
};

/*----------------------------------------------------------------------------*/

/*! \brief Replacement for terminal/terminal.h's color function that will 
 *         ensure the global singleton terminal object has been constructed 
 *         prior to the use of color strings. */
template <typename T> std::string ccolor(T const value) {
  if (Term::DracoTerminal::is_initialized())
    return "\033[" + std::to_string(static_cast<int>(value)) + "m";
  else
    return std::string();
}

} // namespace Term

#endif // dsxx_dracoterminal_hh

//---------------------------------------------------------------------------//
// end of ds++/DracoTerminal.hh
//---------------------------------------------------------------------------//
