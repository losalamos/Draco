//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/OdfmgOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:06:13 2001
 * \brief  OdfmgOpacity class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "OdfmgOpacity.hh"

namespace rtt_cdi
{

/*!
 * \brief Default Opacity() destructor.
 *
 * This is required to correctly release memory when any
 * object derived from OdfmgOpacity is destroyed.
 */
OdfmgOpacity::~OdfmgOpacity()
{
  // empty
}

} // end namespace rtt_cdi


//---------------------------------------------------------------------------//
// end of OdfmgOpacity.cc
//---------------------------------------------------------------------------//
