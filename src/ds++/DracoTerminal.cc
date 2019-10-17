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

// Initialize pointer to zero so that it can be initialized in first call to
// getInstance
Term::DracoTerminal *Term::DracoTerminal::instance = nullptr;

//---------------------------------------------------------------------------//
// end of ds++/DracoTerminal.cc
//---------------------------------------------------------------------------//
