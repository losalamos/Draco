//----------------------------------*-C++-*----------------------------------//
// P1Diffusion.t.cc
// Randy M. Roberts
// Tue Sep 22 15:48:54 1998
//---------------------------------------------------------------------------//
// @> Source for General solution of cell-centered P1 equation.
//---------------------------------------------------------------------------//

#include "P1Diffusion.hh"
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
				 const SP<MS> &spsolver_,
                                 const FieldConstructor &fCtor_)
     : spm(spm_), spsolver(spsolver_), fCtor(fCtor_)
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

    // Get the areas of each face.
    
     fcdsf areas(fCtor);
     spm->get_face_areas(areas);

    // Store the off-diagonal elements of the matrix in a face-centered-
    // discontinous scalar field.
    
     SP<fcdsf> spAOffDiagonal(new fcdsf(fCtor));
     *spAOffDiagonal = -1.0 * areas * DEffOverDeltaL;
    
    // Get the cell volumes.
    
     ccsf volumes(fCtor);
     spm->get_cell_volumes(volumes);

    // Store the diagonal elements of the matrix in a cell-centered
    // scalar field.
    // One part of the diagonal elements include the removal (sigma)
    // mulitplied by the cell volumes.
    // The other part of the diagonal elements includes the negative
    // sum over faces of the off-diagonal elements .

     SP<ccsf> spADiagonal(new ccsf(fCtor));
    
     *spADiagonal = volumes*sigma;
     MT::scatter(*spADiagonal, *spAOffDiagonal, MT::OpSubAssign());

     // Calculate the effective Fprime on the faces due to flux
     // continuity, and the elimination of the face phi unknown.
     // (for internal and boundary faces)
    
     fcdsf FprimeEff(fCtor);
     getFprimeEff(FprimeEff, Fprime, D, alpha, beta, fb);
    
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

#ifdef P13T_MOMENTUM_DEPOSITION

 template<class MT, class MS>
 void
 P1Diffusion<MT,MS>::discFluxToDiscMomentum(DiscMomentumField &result,
					    const DiscFluxField &flux) const
 {
     // This method moves the flux-like field from the DiscFluxField location
     // to the DiscMomentumField location.

     typedef typename MT::fcdvsf NormalsField;

     NormalsField faceNormals(fCtor);
     faceNormals.get_Mesh().get_face_normals(faceNormals);

     // ConnFacesAroundVertices is a class that defines objects
     // that will iterate through a face-centered field around each vertex,
     // before going onto the next vertex's faces.

     typedef typename MT::ConnFacesAroundVertices<NormalsField> ConnNormals;
     typedef typename MT::ConnFacesAroundVertices<const DiscFluxField> ConnFlux;

     const ConnNormals connNormals(faceNormals);
     const ConnFlux connFlux(flux);

     // Make sure that everyone is the same size (in vertices).

     Assert(std::distance(connNormals.begin(), connNormals.end()) ==
	    std::distance(connFlux.begin(), connFlux.end()));
     Assert(std::distance(connNormals.begin(), connNormals.end()) ==
	    std::distance(result.begin(), result.end()));
     
     // Loop over vertices
     
     ConnNormals::const_iterator itNormVertex = connNormals.begin();
     ConnFlux::const_iterator itFluxVertex = connFlux.begin();
     DiscMomentumField::iterator itResult = result.begin();

     while (itResult != result.end())
     {
	 // Make sure that everyone is the same size (in faces per vertex).
	 
	 Assert(std::distance((*itNormVertex).begin(),
			      (*itNormVertex).end()) ==
		std::distance((*itFluxVertex).begin(),
			      (*itFluxVertex).end()));
	 
	 // Loop over faces per vertex
	 
	 ConnNormals::value_type::const_iterator itNormFace =
	     (*itNormVertex).begin();
	 ConnFlux::value_type::const_iterator itFluxFace =
	     (*itFluxVertex).begin();

	 while (itNormFace != (*itNormVertex).end())
	 {
	     typedef NormalsField::value_type NormVector;
	     typedef DiscMomentumField::value_type ResultVector;
	     typedef DiscFluxField::value_type FluxType;

	     ResultVector &res = *itResult;
	     const NormVector &norm = *itNormFace;
	     const FluxType &flux = *itFluxFace;

	     res[0] = norm[0]*flux;
	     res[1] = norm[1]*flux;
	     res[2] = norm[2]*flux;
	     
	     itNormFace++;
	     itFluxFace++;
	 }

	 itResult++;
	 itNormVertex++;
	 itFluxVertex++;
     }
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::dotProduct(DiscKineticEnergyField &result,
                                     const DiscMomentumField &vec1,
                                     const DiscMomentumField &vec2) const
 {
     Assert(result.size() == vec1.size());
     
     DiscMomentumField::const_iterator iv1 = vec1.begin();
     DiscMomentumField::const_iterator iv2 = vec2.begin();
     DiscKineticEnergyField::iterator ir = result.begin();

     while (ir != result.end())
     {
	 typedef DiscMomentumField::value_type vector;	 
	 typedef rtt_traits::vector_traits<vector> vtraits;

	 *ir = vtraits::dot(*iv1, *iv2);

	 ir++;
	 iv1++;
	 iv2++;
     }
 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::dotProduct(DiscKineticEnergyField &KEnergy,
                                     const DiscFluxField &sigmaF,
                                     const DiscMomentumField &velocity) const
 {
     Assert(KEnergy.size() == velocity.size());

     // Move the flux-like field from the DiscFluxField location
     // to the DiscMomentumField location.
     
     DiscMomentumField sigmaFAtMomentum(fCtor);
     discFluxToDiscMomentum(sigmaFAtMomentum, sigmaF);

     // KEnergy is the
     // dot_product((sigma * Flux), velocity)
     // at each vertex.

     dotProduct(KEnergy, sigmaFAtMomentum, velocity);
     
 }

#endif

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::cacheSwappedValues(const fcdsf &D,
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

     bssf DBndry(fCtor);
     bssf deltaLBndry(fCtor);

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

     bssf DBndry(fCtor);
     bssf FprimeBndry(fCtor);
     bssf deltaLBndry(fCtor);

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
    
     SP<Matrix> spMatrix(new Matrix(fCtor, spADiagonal, spAOffDiagonal));

     // Solve the "matrix*phi = b" equations.
    
     spsolver->solve(phi, spMatrix, brhs);

 }

 template<class MT, class MS>
 void P1Diffusion<MT,MS>::getNewFlux(fcdsf &F, const fcdsf &DEffOverDeltaL,
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

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                              end of P1Diffusion.t.cc
//---------------------------------------------------------------------------//
