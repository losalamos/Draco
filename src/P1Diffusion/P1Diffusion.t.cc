//----------------------------------*-C++-*----------------------------------//
// P1Diffusion.t.cc
// Randy M. Roberts
// Tue Sep 22 15:48:54 1998
//---------------------------------------------------------------------------//
// @> Source for General solution of cell-centered P1 equation.
//---------------------------------------------------------------------------//

#include "P1Diffusion/P1Diffusion.hh"
#include "ds++/SP.hh"

// Note:
// Inside class scope MT has been typedef'ed to MeshType,
// and MS has been typedef'ed to MatrixSolver.

namespace rtt_P1Diffusion
{
 // Since rtt_P1Diffusion is a limited namespace
 // I risk putting in this using declaration.
 
 using dsxx::SP;
 
 template<class MT, class MS>
 P1Diffusion<MT,MS>::P1Diffusion(const Diffusion_DB &diffdb,
				 const SP<MT>& spm_,
				 const SP<MS> &spsolver_)
     : spm(spm_), spsolver(spsolver_)
 {
     // empty
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::solve(ccsf &phi, fcdsf &F,
				const fcdsf &D, const ccsf &sigma,
				const ccsf &Q, const fcdsf &Fprime,
				const bssf &alpha, const bssf &beta,
				const bssf &fb) const
 {
     // Get the cell lengths perpendicular to the face.

     fcdsf deltaL(spm);
     spm->get_face_lengths(deltaL);
     
    // Cache swapped values of
    //   D, Fprime, and deltaL

     cacheSwappedValues(D, Fprime, deltaL);
    
     // Calculate the effective D/delta l on the faces due to flux
     // continuity, and the elimination of the face phi unknown.
     // (for internal and boundary faces)
    
     fcdsf DEffOverDeltaL(spm);
     getDEffOverDeltaL(DEffOverDeltaL, D, alpha, beta);

    // Get the areas of each face.
    
     fcdsf areas(spm);
     spm->get_face_areas(areas);

    // Store the off-diagonal elements of the matrix in a face-centered-
    // discontinous scalar field.
    
     SP<fcdsf> spAOffDiagonal = new fcdsf(spm);
     *spAOffDiagonal = -1.0 * areas * DEffOverDeltaL;
    
    // Get the cell volumes.
    
     ccsf volumes(spm);
     spm->get_cell_volumes(volumes);

    // Store the diagonal elements of the matrix in a cell-centered
    // scalar field.
    // One part of the diagonal elements include the removal (sigma)
    // mulitplied by the cell volumes.
    // The other part of the diagonal elements includes the negative
    // sum over faces of the off-diagonal elements .

     SP<ccsf> spADiagonal = new ccsf(spm);
    
     *spADiagonal = volumes*sigma;
     MT::scatter(*spADiagonal, *spAOffDiagonal, MT::OpSubAssign());

     // Calculate the effective Fprime on the faces due to flux
     // continuity, and the elimination of the face phi unknown.
     // (for internal and boundary faces)
    
     fcdsf FprimeEff(spm);
     getFprimeEff(FprimeEff, Fprime, D, alpha, beta, fb);
    
    // The right hand side of the matrix equation is a cell-centered
    // scalar field.

     ccsf brhs(spm);

     // This bracket is for scoping.
     {
	 // Need a temporary face-centered field for continuing the rhs
	 // calculation.
	
	 fcdsf bf(spm);

	 // The right hand side consists the volume multiplied by the source...
	 // The rhs also consists of the negative sum over faces
	 // of the product of effective flux sources and areas.

	 brhs = volumes*Q;
	 bf = areas * FprimeEff;
	 MT::scatter(brhs, bf, MT::OpSubAssign());
     }

     // We no longer need the cached values.
    
     decacheSwappedValues();
    
     // Solve the "matrix*phi = b" equations,
     // where the matrix is comprised of the diagonal and off-diagonal
     // parts.

     solveMatrixEquation(phi, spADiagonal, spAOffDiagonal, brhs);

     // Calculate the new flux using the new values of phi.
    
     getNewFlux(F, DEffOverDeltaL, FprimeEff, phi);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::cacheSwappedValues(const fcdsf &D,
					     const fcdsf &Fprime,
					     const fcdsf &deltaL) const
 {
     spDSwap = new fcdsf(spm);
     MT::swap(*spDSwap, D);

     spFprimeSwap = new fcdsf(spm);
     MT::swap(*spFprimeSwap, Fprime);

     spDeltaLSwap = new fcdsf(spm);
     MT::swap(*spDeltaLSwap, deltaL);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::decacheSwappedValues() const
 {
     spDSwap = SP<fcdsf>();
     spFprimeSwap = SP<fcdsf>();
     spDeltaLSwap = SP<fcdsf>();
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getDEffOverDeltaL(fcdsf &DEffOverDeltaL,
					    const fcdsf &D,
					    const bssf &alpha,
					    const bssf &beta) const
 {
     // Get the cell lengths perpendicular to the face.

     fcdsf deltaL(spm);
     spm->get_face_lengths(deltaL);
    
    // Calculate the effective D/delta l on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (internal faces only)
    
     getDEffOverDeltaLInternal(DEffOverDeltaL, D, deltaL);

     // Calculate the effective D/delta l on the faces due to flux
     // continuity, and the elimination of the face phi unknown.
     // (boundary faces only)
    
     bssf DEffOverDeltaLBndry(spm);
     getDEffOverDeltaLBndry(DEffOverDeltaLBndry, D, deltaL, alpha, beta);

    // Combine the boundary term with interior term.

     MT::gather(DEffOverDeltaL, DEffOverDeltaLBndry, MT::OpAssign());
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getDEffOverDeltaLInternal(fcdsf &DEffOverDeltaL,
						    const fcdsf &D,
						    const fcdsf &deltaL) const
 {
     // Need to calculate the D's and deltaL's from the opposite face.

     const fcdsf &DSwap = getDSwap();
     const fcdsf &deltaLSwap = getDeltaLSwap();
    
     DEffOverDeltaL = 2.0*D*DSwap / (D*deltaLSwap + DSwap*deltaL);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getDEffOverDeltaLBndry(bssf &DEffOverDeltaLBndry,
						 const fcdsf &D,
						 const fcdsf &deltaL,
						 const bssf &alpha,
						 const bssf &beta) const
 {
     // Need to strip off the boundary part of the D's and deltaL's.

     bssf DBndry(spm);
     bssf deltaLBndry(spm);

     MT::gather(DBndry, D, MT::OpAssign());
     MT::gather(deltaLBndry, deltaL, MT::OpAssign());

     DEffOverDeltaLBndry = 2.0*alpha*DBndry /
	 (alpha*deltaLBndry - 2.0*beta*DBndry);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getFprimeEff(fcdsf &FprimeEff,
				       const fcdsf &Fprime,
				       const fcdsf &D,
				       const bssf &alpha,
				       const bssf &beta,
				       const bssf &fb) const
 {
     // Get the cell lengths perpendicular to the face.

     fcdsf deltaL(spm);
     spm->get_face_lengths(deltaL);
    
    // Calculate the effective Fprime on the faces due to flux
    // continuity, and the elimination of the face phi unknown.
    // (internal faces only)
    
     getFprimeEffInternal(FprimeEff, Fprime, D, deltaL);
    
     // Calculate the effective Fprime on the faces due to flux
     // continuity, and the elimination of the face phi unknown.
     // (boundary faces only)
    
     bssf FprimeEffBndry(spm);
     getFprimeEffBndry(FprimeEffBndry, Fprime, D, deltaL, alpha, beta, fb);

    // Combine the boundary term with interior term.

     MT::gather(FprimeEff, FprimeEffBndry, MT::OpAssign());
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getFprimeEffInternal(fcdsf &FprimeEff,
					       const fcdsf &Fprime,
					       const fcdsf &D,
					       const fcdsf &deltaL) const
 {
     // Need to calculate the D's, Fprime's, and deltaL's
     // from the opposite face.

     const fcdsf &DSwap = getDSwap();
     const fcdsf &FprimeSwap = getFprimeSwap();
     const fcdsf &deltaLSwap = getDeltaLSwap();
    
     FprimeEff = (DSwap*deltaL*Fprime - D*deltaLSwap*FprimeSwap) /
	 (D*deltaLSwap + DSwap*deltaL);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getFprimeEffBndry(bssf &FprimeEffBndry,
					    const fcdsf &Fprime,
					    const fcdsf &D,
					    const fcdsf &deltaL,
					    const bssf &alpha,
					    const bssf &beta,
					    const bssf &fb) const
 {
     // Need to strip off the boundary part of the D's, Fprime's, and deltaL's.

     bssf DBndry(spm);
     bssf FprimeBndry(spm);
     bssf deltaLBndry(spm);

     MT::gather(DBndry, D, MT::OpAssign());
     MT::gather(FprimeBndry, Fprime, MT::OpAssign());
     MT::gather(deltaLBndry, deltaL, MT::OpAssign());

     FprimeEffBndry = (deltaLBndry*alpha*FprimeBndry - 2.0*DBndry*fb) /
	 (alpha*deltaLBndry - 2.0*beta*DBndry);
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::
 solveMatrixEquation(ccsf &phi, const SP<const ccsf> &spADiagonal,
		     const SP<const fcdsf> &spAOffDiagonal,
		     const ccsf &brhs) const
 {
     // The sparse matrix can be constructed from the diagonal and
     // off-diagonal elements.
    
     SP<Matrix> spMatrix = new Matrix(spm, spADiagonal, spAOffDiagonal);

     // Solve the "matrix*phi = b" equations.
    
     spsolver->solve(phi, spMatrix, brhs);

 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getNewFlux(fcdsf &F, const fcdsf &DEffOverDeltaL,
				     const fcdsf &FprimeEff,
				     const ccsf &phi) const
 {
     // Get the cell-centered phi onto the faces.

     fcdsf phiFC(spm);
     MT::gather(phiFC, phi, MT::OpAssign());
    
    // Get the values of phi from the opposite side of the face.
    // (Boundary face values will be zero-ed out by the swap.)
    
     fcdsf phiFCSwap(spm);
     MT::swap(phiFCSwap, phiFC);

    // Calculate the new flux.
    
     F = FprimeEff - DEffOverDeltaL*(phiFCSwap - phiFC);
 }

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                              end of P1Diffusion.t.cc
//---------------------------------------------------------------------------//
