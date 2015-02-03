//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/Index_Labeled.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 13:57:58 2005
 * \brief  Header for Index_Labeled.
 * \note   Copyright (C) 2005-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_norms_Index_Labeled_hh
#define rtt_norms_Index_Labeled_hh

#include "ds++/config.h"
#include <string>

#if defined(MSVC)
#   pragma warning (push)
#   pragma warning (disable:4251) // warning C4251: 
#endif

namespace rtt_norms
{

//===========================================================================//
/*!
 * \struct Index_Labeled
 * \brief  A Norms index class that contains processor number and a label.
 */
//===========================================================================//

struct DLL_PUBLIC Index_Labeled 
{
    //! The index.
    size_t index;

    //! The processor number.
    size_t processor;

    //! A string label.
    std::string label;

    // Constructor.  Allow auto-casting from index.
    Index_Labeled(const size_t index_ = 0,
		  const std::string &label_ = "");

    // Equality operator.
    bool operator==(const Index_Labeled &rhs) const;
};

} // end namespace rtt_norms

#if defined(MSVC)
#   pragma warning (pop)
#endif

#endif // rtt_norms_Index_Labeled_hh

//---------------------------------------------------------------------------//
// end of norms/Index_Labeled.hh
//---------------------------------------------------------------------------//
