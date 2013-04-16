//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/cdi_test.hh
 * \author Thomas M. Evans
 * \date   Tue Oct  9 10:51:39 2001
 * \brief  CDI Test help function prototypes.
 * \note   Copyright (C) 2001-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_test_hh__
#define __cdi_test_hh__

#include "ds++/config.h"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_cdi_test
{

//---------------------------------------------------------------------------//
// CHECK COMPUTED VERSUS EXPECTED VALUES
//---------------------------------------------------------------------------//

DLL_PUBLIC bool match(double computedValue, double referenceValue);

DLL_PUBLIC bool match(const std::vector< double > &computedValue, 
	   const std::vector< double > &referenceValue );

DLL_PUBLIC bool match(
    const std::vector< std::vector< double > >& computedValue, 
    const std::vector< std::vector< double > >& referenceValue ); 

DLL_PUBLIC bool match(
    const std::vector< std::vector< std::vector< double > > >& computedValue, 
    const std::vector< std::vector< std::vector< double > > >& referenceValue ); 

} // end namespace rtt_cdi_test

#endif // __cdi_test_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/cdi_test.hh
//---------------------------------------------------------------------------//
