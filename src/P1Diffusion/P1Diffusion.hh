//----------------------------------*-C++-*----------------------------------//
// P1Diffusion.hh
// Randy M. Roberts
// Tue Sep 15 11:08:28 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_P1Diffusion_hh__
#define __P1Diffusion_P1Diffusion_hh__

#include "P1Momentum.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "diffusion/Diffusion_DB.hh"
#include "traits/MatrixFactoryTraits.hh"

// Forward reference
namespace rtt_diffusion
{
template<class MT> class P1Matrix;
}

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

 template<class MT, class MS, bool HASVELOCITY=true>
 class P1Diffusion
 {

     // NESTED CLASSES AND TYPEDEFS

   private:
     
     typedef MS MatrixSolver;
     typedef typename MS::Matrix Matrix;
     typedef typename MT::fcdsf fcdsf;
     typedef typename MT::ccsf ccsf;
     typedef typename MT::bssf bssf;
     typedef typename MT::FieldConstructor FieldConstructor;

     typedef rtt_traits::MatrixFactoryTraits<Matrix> MatFacTraits;
     typedef rtt_diffusion::P1Matrix<MT> P1Matrix;
     
   public:

     typedef MT MeshType;
     typedef P1Momentum<MT,HASVELOCITY> P1Momentum;
     typedef ccsf  IntensityField;
     typedef fcdsf FluxField;
     typedef fcdsf DiscFluxField;
     typedef fcdsf DiffCoefField;
     typedef typename P1Momentum::MomentumField MomentumField;
     typedef typename P1Momentum::DiscMomentumField DiscMomentumField;
     typedef typename P1Momentum::MomentumComponentField
                                                    MomentumComponentField;
     typedef typename P1Momentum::DiscMomentumComponentField
                                                    DiscMomentumComponentField;
     typedef typename P1Momentum::DiscKineticEnergyField
                                                    DiscKineticEnergyField;

     // DATA
    
     rtt_dsxx::SP<MeshType> spm;
     rtt_dsxx::SP<MatrixSolver> spsolver;
     rtt_dsxx::SP<P1Momentum> spmomentum;
     FieldConstructor fCtor;
     typename MatFacTraits::PreComputedState preComputedMatrixState;
     bool jacobiScale;
     
     // Cache the swapped values to avoid too much communication.
     
     mutable rtt_dsxx::SP<fcdsf> spDSwap;
     mutable rtt_dsxx::SP<fcdsf> spFprimeSwap;
     mutable rtt_dsxx::SP<fcdsf> spDeltaLSwap;

   public:

     // CREATORS
    
     P1Diffusion(const rtt_diffusion::Diffusion_DB &diffdb,
		 const rtt_dsxx::SP<MeshType>& spm_,
		 const rtt_dsxx::SP<MatrixSolver> &spsolver_,
                 const FieldConstructor &fCtor_);

     // MANIPULATORS
    
     // ACCESSORS

     void solve(ccsf &phi, fcdsf &F, const fcdsf &D, const ccsf &sigma,
		const ccsf &Q, const fcdsf &Fprime, const bssf &alpha,
		const bssf &beta, const bssf &fb) const;

     typename ccsf::value_type integrateOverVolume(const ccsf &field) const;

     typename fcdsf::value_type integrateOverBoundary(const fcdsf &field) const;

    void discFluxToDiscMomentum(DiscMomentumField &result,
				const DiscFluxField &flux) const
    {
	spmomentum->discFluxToDiscMomentum(result, flux);
    }

    void dotProduct(DiscKineticEnergyField &result,
		    const DiscMomentumField &vec1,
		    const DiscMomentumField &vec2) const
    {
	spmomentum->dotProduct(result, vec1, vec2);
    }
     
    void dotProduct(DiscKineticEnergyField &KEnergy,
		    const DiscFluxField &sigmaF,
		    const DiscMomentumField &velocity) const
    {
	spmomentum->dotProduct(KEnergy, sigmaF, velocity);
    }

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

     void solveMatrixEquation(ccsf &phi, P1Matrix &p1Mat,
			      const ccsf &brhs) const;
     
     void getNewFlux(fcdsf &F, const fcdsf &DEffOverDeltaL,
		     const fcdsf &FprimeEff,
		     const ccsf &phi) const;

     P1Matrix getP1Matrix(const fcdsf &D, const fcdsf &DEffOverDeltaL,
			  const ccsf &sigma, const bssf &alpha,
			  const bssf &beta) const;

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
