//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.i.hh
 * \author Kelly Thompson
 * \date   Mon Nov  8 14:23:03 2004
 * \brief  Member definitions of class QuadServices
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_QuadServices_i_hh
#define quadrature_QuadServices_i_hh

#include <iostream>
#include "QuadServices.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief kronecker_delta
 * 
 * Return 1 if test_value == offset, otherwise return 0;
 * 
 * \param test_value
 * \param offset
 * \return 1 if test_value == offset, otherwise return 0;
 *
 * \bug move to ds++
 */
template< typename T >
T QuadServices::kronecker_delta( T const test_value,
				 T const offset ) const
{ 
    return (test_value==offset) ? 1 : 0; 
}

//---------------------------------------------------------------------------//
/*! 
 * \brief factorial
 * 
 * This is a recursive function that computes k!
 *
 * \param k 
 * \return k!
 *
 * \bug move to ds++
 */
template< typename T >
T QuadServices::factorial( T const k ) const
{
    if( k <= 1 ) 
	return 1;
    else 
	return  k * factorial(k-1);
}

//---------------------------------------------------------------------------//

template< typename T >
void QuadServices::print_matrix( std::string    const & matrix_name,
				 std::vector<T> const & x,
				 std::vector<unsigned> const & dims ) const
{
    using std::cout;
    using std::endl;
    using std::string;

    unsigned pad_len( matrix_name.length()+2 );
    string padding( pad_len, ' ' );
    cout << matrix_name << " =";
    for( unsigned i=0; i<dims[0]; ++i )
    {
	if( i != 0 ) cout << padding;
	cout << "{ ";
	for( unsigned j=0; j<dims[1]-1; ++j )
	    cout << x[j+dims[0]*i] << ", ";
	cout << x[dims[1]-1+dims[0]*i] << " }." << endl;
    }
    return;
}


} // end namespace rtt_quadrature

#endif // quadrature_QuadServices_i_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.i.hh
//---------------------------------------------------------------------------//
