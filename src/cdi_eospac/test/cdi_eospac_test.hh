//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/cdi_eospac_test.hh
 * \author Kelly Thompson
 * \date   Mon Apr 2 14:15:57 2001
 * \brief  Header file for cdi_eospac_test.cc and tEospac.cc
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_test_hh__
#define __cdi_eospac_test_hh__

#include <vector>
#include <string>

namespace rtt_cdi_eospac_test
{

//---------------------------------------------------------------------------//
// DATA EQUIVALENCE FUNCTIONS USED FOR TESTING
//---------------------------------------------------------------------------//

bool match( const double computedValue, const double referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< double >& computedValue, 
	   const std::vector< double >& referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< std::vector<double> >& computedValue, 
	   const std::vector< std::vector<double> >& referenceValue );

} // end namespace rtt_cdi_eospac_test

#endif // __cdi_eospac_test_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/test/tEospac.hh
//---------------------------------------------------------------------------//

