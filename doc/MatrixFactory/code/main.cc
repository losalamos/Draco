#include "JoubertMat.hh"
#include "JoubertMatTraits.hh"
#include "DenseMatrixRep.hh"

#include "ds++/SP.hh"

#include <iostream>

using namespace rtt_MatrixFactory;

int main()
{
    const double matdata[] = {0.0, 1.0, 2.0, 3.0,
			      0.0, 5.0, 0.0, 7.0,
			      0.0, 0.0, 0.0, 11.0};
    const int matdataSize = sizeof(matdata)/sizeof(double);

    using JoubertMat::JoubertMat;
    
    dsxx::SP<JoubertMat> spJoubertMat;
    JoubertMat *pJoubertMat;

    // Scoping braces
    {

	// Create a Dense matrix to test out the MatrixFactoryTraits
	
	DenseMatrixRep denseMat(3, 4, matdata+0, matdata+matdataSize);
	std::cerr << denseMat << std::endl;

	// Use the factory traits to create a smart pointer to a Joubert matrix.
	
	spJoubertMat = MatrixFactoryTraits<JoubertMat>::create(denseMat);

	// Use the factory traits to create a dumb pointer to a Joubert matrix.
	
	pJoubertMat = MatrixFactoryTraits<JoubertMat>::create(denseMat);
    }

    // Print out the Joubert matrices.
    
    std::cerr << *spJoubertMat << std::endl;
    std::cerr << *pJoubertMat << std::endl;

    // Delete the Joubert matrix via the dumb pointer.
    
    std::cerr << "Deleting pJoubertMat" << std::endl;
    delete pJoubertMat;

    // Delete the Joubert matrix via the smart pointer.
    
    std::cerr << "Zeroing spJoubertMat" << std::endl;
    spJoubertMat = 0;
    
    std::cerr << "All's well." << std::endl;

    return 0;
}
