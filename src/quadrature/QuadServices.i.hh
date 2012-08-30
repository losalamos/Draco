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
#include <iomanip>

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Pretty print vector<T> as 2D matrix.
 * \param matrix_name A string that will be printed as an identifier for this
 *        matrix.
 * \param x The vector<T> that we want to print in a 2D format.
 * \param dims A length 2 vector that contains the dimensions of x.
 * \return void
 */
template< typename T >
void QuadServices::print_matrix( std::string    const & matrix_name,
				 std::vector<T> const & x,
				 std::vector<unsigned> const & dims ) const
{
    using std::cout;
    using std::endl;
    using std::string;

    Require( dims[0]*dims[1] == x.size() );

    unsigned pad_len( matrix_name.length()+2 );
    string padding( pad_len, ' ' );
    cout << matrix_name << " =";
    // row
    for( unsigned i=0; i<dims[1]; ++i )
    {
	if( i != 0 ) cout << padding;

	cout << "{ ";

	for( unsigned j=0; j<dims[0]-1; ++j )
	    cout << setprecision(10) << x[j+dims[0]*i] << ", ";

	cout << setprecision(10) << x[dims[0]-1+dims[0]*i] << " }." << endl;
    }
    cout << endl;
    return;
}


} // end namespace rtt_quadrature

#endif // quadrature_QuadServices_i_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.i.hh
//---------------------------------------------------------------------------//
