/*-----------------------------------*-C++-*---------------------------------*/
/* main.cc */
/* Randy M. Roberts */
/* Thu May 27 16:12:42 1999 */
/*---------------------------------------------------------------------------*/
/* @> Main program to demonstrate MatrixFactoryTraits mechanism */
/*---------------------------------------------------------------------------*/

#include "JoubertMat.hh"
#include "JoubertMatTraits.hh"
#include "DenseMatrixRep.hh"

#include "ds++/SP.hh"

#include <iostream>

using namespace rtt_MatrixFactory;

// use unnamed namespace to keep from polluting the global namespace

namespace
{

// The function doit is templated on an unknown Matrix type.
// This shows the generic use of the MatrixFactoryTraits class.

template<class Matrix>
void doit()
{
    const double matdata[] = {0.0, 1.0, 2.0, 3.0,
			      0.0, 5.0, 0.0, 7.0,
			      0.0, 0.0, 0.0, 11.0};
    const int matdataSize = sizeof(matdata)/sizeof(double);

    dsxx::SP<Matrix> spMatrix;
    Matrix *pMatrix;

    // Scoping braces to destroy denseMat after no longer needed.
    {

	// Create a Dense matrix to test out the MatrixFactoryTraits
	
	DenseMatrixRep denseMat(3, 4, matdata+0, matdata+matdataSize);
	std::cerr << denseMat << std::endl;

	// Use the factory traits to create a smart pointer to an
	// unknown matrix.
	
	spMatrix = MatrixFactoryTraits<Matrix>::create(denseMat);

	// Use the factory traits to create a dumb pointer to an
	// unknown matrix.
	
	pMatrix = MatrixFactoryTraits<Matrix>::create(denseMat);
    }

    // Print out the unknown matrices.
    
    std::cerr << *spMatrix << std::endl;
    std::cerr << *pMatrix << std::endl;

    // Delete the unknown matrix via the dumb pointer.
    
    std::cerr << "Deleting pMatrix" << std::endl;
    delete pMatrix;

    // Delete the unknown matrix via the smart pointer.
    
    std::cerr << "Zeroing spMatrix" << std::endl;
    spMatrix = 0;
    
    std::cerr << "All's well." << std::endl;

}

} // end namespace

int main()
{
    using JoubertNS::JoubertMat;

    // Call doit for a Joubert Matrix.
    
    doit<JoubertMat>();

    return 0;
}

/*---------------------------------------------------------------------------*/
/*    end of main.cc */
/*---------------------------------------------------------------------------*/
