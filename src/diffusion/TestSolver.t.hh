//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diffusion/TestSolver.t.hh
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:26:24 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestSolver.hh"
#include "P1Matrix.hh"
#include "ds++/SP.hh"

namespace rtt_diffusion
{

template<class MT, class Solver>
void TestSolver<MT,Solver>::runTest()
{
    setPassed(true);

    os() << "---------------------------" << std::endl;
    os() << "Without Row/Column Scaling." << std::endl;
    os() << "---------------------------" << std::endl;
    runTest(false);

    os() << "---------------------------" << std::endl;
    os() << "With Row/Column Scaling." << std::endl;
    os() << "---------------------------" << std::endl;
    runTest(true);
}

template<class MT, class Solver>
void TestSolver<MT,Solver>::runTest(bool jacobiScale)
{
    using rtt_dsxx::SP;
    
    if (jacobiScale)
	diff_db.pc_meth = 1;
    else
	diff_db.pc_meth = 0;

    P1Matrix<MT> p1Mat = createP1Matrix();

    MatFacTraits::PreComputedState preComputedMatrixState =
	MatFacTraits::preComputeState(fCtor, p1Mat.diagonal().get_Mesh());
	
    // Create the solver's matrix with a matrix factor traits method.
     
    SP<Matrix> spMatrix(MatFacTraits::create(p1Mat, preComputedMatrixState));

    ccsf phi(fCtor);
    ccsf brhs(fCtor);
    ccsf phi0(fCtor);

    for (int i=0; i<brhs.size(); i++)
	phi0[i] = i;

    multiply(brhs, *spMatrix, phi0);
    // spMatrix->multiply(brhs, phi0);
    
    ccsf D1_2(fCtor);
    if (jacobiScale)
    {
	D1_2 = p1Mat.diagonal();
	D1_2 = sqrt(fabs(D1_2));
	
	p1Mat.jacobiScale();

	// Create the new scaled matrix
	spMatrix = MatFacTraits::create(p1Mat, preComputedMatrixState);
    
	brhs = brhs / D1_2;
	solver.solve(phi, spMatrix, brhs);
	phi = phi / D1_2;
    }
    else
    {
	solver.solve(phi, spMatrix, brhs);
    }

    double error = 0.0;
    double l2nPhi0 = 0.0;
    for (int i=0; i<phi.size(); i++)
    {
	error += (phi[i]-phi0[i])*(phi[i]-phi0[i]);
	l2nPhi0 += phi0[i]*phi0[i];
    }

    C4::gsum(error);
    C4::gsum(l2nPhi0);
     
    error /= l2nPhi0;
    error = std::sqrt(error);
     
    os() << "error: " << error << std::endl;
    os() << "tolerance: " << tolerance << std::endl;

    testAssert(error < tolerance, "error < tolerance", __FILE__, __LINE__);
    
    // os() << std::endl << phi << std::endl;
}    

template<class MT, class Solver>
P1Matrix<MT> TestSolver<MT,Solver>::createP1Matrix() const
{
    ccsf diag(fCtor);
    
    double d = 0;
    for (ccsf::iterator dit = diag.begin(); dit != diag.end(); ++dit)
    {
	*dit = -6 + d/diag.size(); ++d;
    }
    
    fcdsf offDiag(fCtor);
    offDiag = 1.0;

    // The sparse matrix can be constructed from the diagonal and
    // off-diagonal elements.

    return P1Matrix<MT>(fCtor, diag, offDiag);
}

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                        end of diffusion/TestSolver.t.hh
//---------------------------------------------------------------------------//
