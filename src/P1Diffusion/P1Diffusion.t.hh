//----------------------------------*-C++-*----------------------------------//
// P1Diffusion.t.cc
// Randy M. Roberts
// Tue Sep 22 15:48:54 1998
//---------------------------------------------------------------------------//
// @> Source for General solution of cell-centered P1 equation.
//---------------------------------------------------------------------------//

#include "P1Diffusion.hh"
#include "ds++/SP.hh"
#include "diffusion/P1Matrix.hh"
#include "traits/MatrixFactoryTraits.hh"
#include <limits>
#include <cmath>

// Note:
// Inside class scope MT has been typedef'ed to MeshType,
// and MS has been typedef'ed to MatrixSolver.

namespace rtt_P1Diffusion
{
// Since rtt_P1Diffusion is a limited namespace
// I risk putting in this using declaration.
 
using rtt_dsxx::SP;
 
template<class MT, class MS, bool HV>
P1Diffusion<MT,MS,HV>::P1Diffusion(const rtt_diffusion::Diffusion_DB &diffdb,
				   const SP<MT>& spm_,
				   const SP<MS> &spsolver_,
				   const FieldConstructor &fCtor_)
    : spm(spm_), spsolver(spsolver_), fCtor(fCtor_),
      preComputedMatrixState(MatFacTraits::preComputeState(fCtor_, *spm_)),
      spmomentum(new P1Momentum(fCtor_)), jacobiScale(diffdb.pc_meth == 1)
{
    // empty
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::solve(ccsf &phi, fcdsf &F,
				  const fcdsf &D, const ccsf &sigma,
				  const ccsf &Q, const fcdsf &Fprime,
				  const bssf &alpha, const bssf &beta,
				  const bssf &fb)
{
    assembleTimer_m.start();
    
    // Get the cell lengths perpendicular to the face.

    fcdsf deltaL(fCtor);
    spm->get_face_lengths(deltaL);
     
    // Cache swapped values of
    //   D, Fprime, and deltaL

    cacheSwappedValues(D, Fprime, deltaL);

    // Calculate the effective D/delta l on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (for internal and boundary faces)
    
    fcdsf DEffOverDeltaL(fCtor);
    getDEffOverDeltaL(DEffOverDeltaL, D, alpha, beta);

    // Get the Matrix

    P1Matrix p1Mat = getP1Matrix(D, DEffOverDeltaL, sigma, alpha, beta);    

    // Calculate the effective Fprime on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (for internal and boundary faces)
    
    fcdsf FprimeEff(fCtor);
    getFprimeEff(FprimeEff, Fprime, D, alpha, beta, fb);

    assembleTimer_m.stop();

    rhsTimer_m.start();
    
    // The right hand side of the matrix equation is a cell-centered
    // scalar field.

    ccsf brhs(fCtor);

    // This bracket is for scoping.
    {
	// Need a temporary face-centered field for continuing the rhs
	// calculation.
	
	fcdsf bf(fCtor);

	// The right hand side consists the volume multiplied by the source...
	// The rhs also consists of the negative sum over faces
	// of the product of effective flux sources and areas.

	// Get the cell volumes.
    
	ccsf volumes(fCtor);
	spm->get_cell_volumes(volumes);

	brhs = volumes*Q;

	// Get the areas of each face.
    
	fcdsf areas(fCtor);
	spm->get_face_areas(areas);
	
	bf = areas * FprimeEff;
	MT::scatter(brhs, bf, MT::OpSubAssign());
    }

    rhsTimer_m.stop();
    
    // We no longer need the cached values.
    
    decacheSwappedValues();

    // Solve the "matrix*phi = b" equations,
    // where the matrix is comprised of the diagonal and off-diagonal
    // parts.

    matSolverTimer_m.start();
    solveMatrixEquation(phi, p1Mat, brhs);
    matSolverTimer_m.stop();

    // Calculate the new flux using the new values of phi.

    newFluxTimer_m.start();
    getNewFlux(F, DEffOverDeltaL, FprimeEff, phi);
    newFluxTimer_m.stop();
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::cacheSwappedValues(const fcdsf &D,
					       const fcdsf &Fprime,
					       const fcdsf &deltaL) const
{
    spDSwap = new fcdsf(fCtor);
    MT::swap_faces(*spDSwap, D);

    spFprimeSwap = new fcdsf(fCtor);
    MT::swap_faces(*spFprimeSwap, Fprime);

    spDeltaLSwap = new fcdsf(fCtor);
    MT::swap_faces(*spDeltaLSwap, deltaL);
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::decacheSwappedValues() const
{
    spDSwap = SP<fcdsf>();
    spFprimeSwap = SP<fcdsf>();
    spDeltaLSwap = SP<fcdsf>();
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getDEffOverDeltaL(fcdsf &DEffOverDeltaL,
					      const fcdsf &D,
					      const bssf &alpha,
					      const bssf &beta) const
{
    // Get the cell lengths perpendicular to the face.

    fcdsf deltaL(fCtor);
    spm->get_face_lengths(deltaL);
    
    // Calculate the effective D/delta l on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (internal faces only)
    
    getDEffOverDeltaLInternal(DEffOverDeltaL, D, deltaL);

    // Calculate the effective D/delta l on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (boundary faces only)
    
    bssf DEffOverDeltaLBndry(fCtor);
    getDEffOverDeltaLBndry(DEffOverDeltaLBndry, D, deltaL, alpha, beta);

    // Combine the boundary term with interior term.

    MT::gather(DEffOverDeltaL, DEffOverDeltaLBndry, MT::OpAssign());
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getDEffOverDeltaLInternal(fcdsf &DEffOverDeltaL,
						      const fcdsf &D,
						      const fcdsf &deltaL) const
{
    // Need to calculate the D's and deltaL's from the opposite face.

    const fcdsf &DSwap = getDSwap();
    const fcdsf &deltaLSwap = getDeltaLSwap();

    static const typename fcdsf::value_type
	floor = std::numeric_limits<typename fcdsf::value_type>::min();
    
    DEffOverDeltaL = 2.0*D*DSwap / max((D*deltaLSwap + DSwap*deltaL), floor);
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getDEffOverDeltaLBndry(bssf &DEffOverDeltaLBndry,
						   const fcdsf &D,
						   const fcdsf &deltaL,
						   const bssf &alpha,
						   const bssf &beta) const
{
    // Need to strip off the boundary part of the D's and deltaL's.

    bssf DBndry(fCtor);
    bssf deltaLBndry(fCtor);

    MT::gather(DBndry, D, MT::OpAssign());
    MT::gather(deltaLBndry, deltaL, MT::OpAssign());

    static const typename fcdsf::value_type
	floor = std::numeric_limits<typename fcdsf::value_type>::min();

    bssf denom(fCtor);
    denom = alpha*deltaLBndry - 2.0*beta*DBndry;
    for (bssf::iterator dit = denom.begin(); dit != denom.end(); ++dit)
	*dit = std::fabs(*dit) > floor ? *dit : floor;

    DEffOverDeltaLBndry = 2.0*alpha*DBndry / denom;
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getFprimeEff(fcdsf &FprimeEff,
					 const fcdsf &Fprime,
					 const fcdsf &D,
					 const bssf &alpha,
					 const bssf &beta,
					 const bssf &fb) const
{
    // Get the cell lengths perpendicular to the face.

    fcdsf deltaL(fCtor);
    spm->get_face_lengths(deltaL);
    
    // Calculate the effective Fprime on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (internal faces only)
    
    getFprimeEffInternal(FprimeEff, Fprime, D, deltaL);
    
    // Calculate the effective Fprime on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (boundary faces only)
    
    bssf FprimeEffBndry(fCtor);
    getFprimeEffBndry(FprimeEffBndry, Fprime, D, deltaL, alpha, beta, fb);

    // Combine the boundary term with interior term.

    MT::gather(FprimeEff, FprimeEffBndry, MT::OpAssign());
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getFprimeEffInternal(fcdsf &FprimeEff,
						 const fcdsf &Fprime,
						 const fcdsf &D,
						 const fcdsf &deltaL) const
{
    // Need to calculate the D's, Fprime's, and deltaL's
    // from the opposite face.

    const fcdsf &DSwap = getDSwap();
    const fcdsf &FprimeSwap = getFprimeSwap();
    const fcdsf &deltaLSwap = getDeltaLSwap();
    
    static const typename fcdsf::value_type
	floor = std::numeric_limits<typename fcdsf::value_type>::min();

    FprimeEff = (DSwap*deltaL*Fprime - D*deltaLSwap*FprimeSwap) /
	max((D*deltaLSwap + DSwap*deltaL), floor);
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getFprimeEffBndry(bssf &FprimeEffBndry,
					      const fcdsf &Fprime,
					      const fcdsf &D,
					      const fcdsf &deltaL,
					      const bssf &alpha,
					      const bssf &beta,
					      const bssf &fb) const
{
    // Need to strip off the boundary part of the D's, Fprime's, and deltaL's.

    bssf DBndry(fCtor);
    bssf FprimeBndry(fCtor);
    bssf deltaLBndry(fCtor);

    MT::gather(DBndry, D, MT::OpAssign());
    MT::gather(FprimeBndry, Fprime, MT::OpAssign());
    MT::gather(deltaLBndry, deltaL, MT::OpAssign());

    static const typename fcdsf::value_type
	floor = std::numeric_limits<typename fcdsf::value_type>::min();

    bssf denom(fCtor);
    denom = alpha*deltaLBndry - 2.0*beta*DBndry;
    for (bssf::iterator dit = denom.begin(); dit != denom.end(); ++dit)
	*dit = std::fabs(*dit) > floor ? *dit : floor;

    FprimeEffBndry = (deltaLBndry*alpha*FprimeBndry - 2.0*DBndry*fb) /
	denom;
}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::solveMatrixEquation(ccsf &phi, P1Matrix &p1Mat,
						const ccsf &brhs)
{
    ccsf D1_2(fCtor);
    if (jacobiScale)
    {
	D1_2 = p1Mat.diagonal();
	D1_2 = sqrt(fabs(D1_2));
	
	p1Mat.jacobiScale();
    }

    // Create the solver's matrix with a matrix factor traits method.
     
    SP<Matrix> spMatrix(MatFacTraits::create(p1Mat, preComputedMatrixState));

     // Solve the "matrix*phi = b" equations.
    
    if (jacobiScale)
    {
	ccsf brhsTilde(fCtor);
	brhsTilde = brhs / D1_2;

	spsolver->solve(phi, spMatrix, brhsTilde);
	
	phi = phi / D1_2;
    }
    else
    {
	spsolver->solve(phi, spMatrix, brhs);
    }

}

template<class MT, class MS, bool HV>
void P1Diffusion<MT,MS,HV>::getNewFlux(fcdsf &F, const fcdsf &DEffOverDeltaL,
				       const fcdsf &FprimeEff,
				       const ccsf &phi) const
{
    // Get the cell-centered phi onto the faces.

    fcdsf phiFC(fCtor);
    MT::gather(phiFC, phi, MT::OpAssign());
    
    // Get the values of phi from the opposite side of the face.
    // (Boundary face values will be zero-ed out by the swap.)
    
    fcdsf phiFCSwap(fCtor);
    MT::swap_faces(phiFCSwap, phiFC);

    // Calculate the new flux.
    
    F = FprimeEff - DEffOverDeltaL*(phiFCSwap - phiFC);
}

template<class MT, class MS, bool HV>
typename P1Diffusion<MT,MS,HV>::ccsf::value_type
P1Diffusion<MT,MS,HV>::integrateOverVolume(const ccsf &field) const
{
    // Get the cell volumes.
    
    ccsf result_field(fCtor);
    spm->get_cell_volumes(result_field);

    result_field *= field;

    // return sum over all cells.
    
    return MT::sum(result_field);
}

template<class MT, class MS, bool HV>
typename P1Diffusion<MT,MS,HV>::fcdsf::value_type
P1Diffusion<MT,MS,HV>::integrateOverBoundary(const fcdsf &field) const
{
    // Get the cell volumes.
    
    fcdsf result_fcdsf(fCtor);
    spm->get_face_areas(result_fcdsf);

    result_fcdsf *= field;

    // strip out the boundary of area*field into a boundary field
    
    bssf result_bssf(fCtor);
    MT::gather(result_bssf, result_fcdsf, MT::OpAssign());

    // return the sum over all boundary faces.
    
    return MT::sum(result_bssf);
}

template<class MT, class MS, bool HV>
P1Diffusion<MT,MS,HV>::P1Matrix
P1Diffusion<MT,MS,HV>::getP1Matrix(const fcdsf &D, const fcdsf &DEffOverDeltaL,
				   const ccsf &sigma, const bssf &alpha,
				   const bssf &beta) const
{
    // Get the areas of each face.
    
    fcdsf areas(fCtor);
    spm->get_face_areas(areas);

    // Store the off-diagonal elements of the matrix in a face-centered-
    // discontinous scalar field.
    
    fcdsf AOffDiagonal(fCtor);
    AOffDiagonal = -1.0 * areas * DEffOverDeltaL;
    
    // Get the cell volumes.
    
    ccsf volumes(fCtor);
    spm->get_cell_volumes(volumes);

    // Store the diagonal elements of the matrix in a cell-centered
    // scalar field.
    // One part of the diagonal elements include the removal (sigma)
    // mulitplied by the cell volumes.
    // The other part of the diagonal elements includes the negative
    // sum over faces of the off-diagonal elements .

    ccsf ADiagonal(fCtor);
    
    ADiagonal = volumes*sigma;
    MT::scatter(ADiagonal, AOffDiagonal, MT::OpSubAssign());

    // The sparse matrix can be constructed from the diagonal and
    // off-diagonal elements.

    return P1Matrix(fCtor, ADiagonal, AOffDiagonal);
}

template<class MT, class MS, bool HV>
const std::list<std::pair<std::string,rtt_stopwatch::Timer *> >
P1Diffusion<MT,MS,HV>::timers()
{
    typedef std::pair<std::string,rtt_stopwatch::Timer *> pair;
    std::list<pair> list;
    list.push_back(pair("assemble", &assembleTimer_m));
    list.push_back(pair("rhs", &rhsTimer_m));
    list.push_back(pair("mat solver", &matSolverTimer_m));
    list.push_back(pair("new flux", &newFluxTimer_m));
    return list;
}

template<class MT, class MS, bool HV>
const std::list<std::pair<std::string, const rtt_stopwatch::Timer *> >
P1Diffusion<MT,MS,HV>::timers() const
{
    typedef std::pair<std::string,const rtt_stopwatch::Timer *> pair;
    std::list<pair> list;
    list.push_back(pair("assemble", &assembleTimer_m));
    list.push_back(pair("rhs", &rhsTimer_m));
    list.push_back(pair("mat solver", &matSolverTimer_m));
    list.push_back(pair("new flux", &newFluxTimer_m));
    return list;
}

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                              end of P1Diffusion.t.cc
//---------------------------------------------------------------------------//
