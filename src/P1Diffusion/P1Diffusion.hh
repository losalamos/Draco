//----------------------------------*-C++-*----------------------------------//
// P1Diffusion.hh
// Randy M. Roberts
// Tue Sep 15 11:08:28 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_P1Diffusion_hh__
#define __P1Diffusion_P1Diffusion_hh__

#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "diffusion/Diffusion_DB.hh"

namespace rtt_P1Diffusion
{
    
 //===========================================================================//
 // class P1Diffusion - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class MT, class MS>
 class P1Diffusion
 {

     // NESTED CLASSES AND TYPEDEFS

   private:
     
     typedef MS MatrixSolver;
     typedef typename MS::Matrix Matrix;
     typedef typename MT::fcdsf fcdsf;
     typedef typename MT::ccsf ccsf;
     typedef typename MT::bssf bssf;
#ifdef P13T_MOMENTUM_DEPOSITION
     typedef typename MT::ncvsf ncvsf;
     typedef typename MT::vcvsf vcvsf;
     typedef typename MT::vcsf vcsf;
#endif
     typedef typename MT::FieldConstructor FieldConstructor;

   public:

     typedef MT MeshType;
     typedef fcdsf FluxField;
     typedef fcdsf DiscFluxField;
#ifdef P13T_MOMENTUM_DEPOSITION
     typedef ncvsf MomentumField;
     typedef vcvsf DiscMomentumField;
     typedef vcsf DiscKineticEnergyField;
#endif
     typedef fcdsf DiffCoefField;

     // DATA
    
     dsxx::SP<MeshType> spm;
     dsxx::SP<MatrixSolver> spsolver;
     FieldConstructor fCtor;

     // Cache the swapped values to avoid too much communication.
     
     mutable dsxx::SP<fcdsf> spDSwap;
     mutable dsxx::SP<fcdsf> spFprimeSwap;
     mutable dsxx::SP<fcdsf> spDeltaLSwap;

   public:

     // CREATORS
    
     P1Diffusion(const Diffusion_DB &diffdb, const dsxx::SP<MeshType>& spm_,
		 const dsxx::SP<MatrixSolver> &spsolver_,
                 const FieldConstructor &fCtor_);

     // MANIPULATORS
    
     // ACCESSORS

     void solve(ccsf &phi, fcdsf &F, const fcdsf &D, const ccsf &sigma,
		const ccsf &Q, const fcdsf &Fprime, const bssf &alpha,
		const bssf &beta, const bssf &fb) const;

#ifdef P13T_MOMENTUM_DEPOSITION
     void discFluxToDiscMomentum(DiscMomentumField &result,
				 const DiscFluxField &flux) const;

     void dotProduct(DiscKineticEnergyField &result,
		     const DiscMomentumField &vec1,
		     const DiscMomentumField &vec2) const;
     
     void dotProduct(DiscKineticEnergyField &KEnergy,
                     const DiscFluxField &sigmaF,
                     const DiscMomentumField &velocity) const;
#endif

   private:
    
     // IMPLEMENTATION

     void cacheSwappedValues(const fcdsf &D,
			     const fcdsf &Fprime,
			     const fcdsf &deltaL) const;
     
     void decacheSwappedValues() const;
     
     void getDEffOverDeltaL(fcdsf &DEffOverDeltaL,
			    const fcdsf &D,
			    const bssf &alpha,
			    const bssf &beta) const;

     void getDEffOverDeltaLInternal(fcdsf &DEffOverDeltaL,
				    const fcdsf &D,
				    const fcdsf &deltaL) const;
     
     void getDEffOverDeltaLBndry(bssf &DEffOverDeltaLBndry,
				 const fcdsf &D,
				 const fcdsf &deltaL,
				 const bssf &alpha,
				 const bssf &beta) const;
     
     void getFprimeEff(fcdsf &FprimeEff,
		       const fcdsf &Fprime,
		       const fcdsf &D,
		       const bssf &alpha,
		       const bssf &beta,
		       const bssf &fb) const;

     void getFprimeEffInternal(fcdsf &FprimeEff,
			       const fcdsf &Fprime,
			       const fcdsf &D,
			       const fcdsf &deltaL) const;

     void getFprimeEffBndry(bssf &FprimeEffBndry,
			    const fcdsf &Fprime,
			    const fcdsf &D,
			    const fcdsf &deltaL,
			    const bssf &alpha,
			    const bssf &beta,
			    const bssf &fb) const;

     void solveMatrixEquation(ccsf &phi,
			      const dsxx::SP<const ccsf> &spADiagonal,
			      const dsxx::SP<const fcdsf> &spAOffDiagonal,
			      const ccsf &brhs) const;
     
     void getNewFlux(fcdsf &F, const fcdsf &DEffOverDeltaL,
		     const fcdsf &FprimeEff,
		     const ccsf &phi) const;

     const fcdsf &getDSwap() const
     {
	 Assert(spDSwap);
	 return *spDSwap;
     }
     
     const fcdsf &getFprimeSwap() const
     {
	 Assert(spFprimeSwap);
	 return *spFprimeSwap;
     }
     
     const fcdsf &getDeltaLSwap() const {
	 Assert(spDeltaLSwap);
	 return *spDeltaLSwap;
     }
 };

}  // namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_P1Diffusion_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/P1Diffusion.hh
//---------------------------------------------------------------------------//
