//----------------------------------*-C++-*----------------------------------//
// P13T.hh
// Randy M. Roberts
// Wed Mar 11 11:18:52 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_P13T_hh__
#define __3T_P13T_hh__

#include "ds++/SP.hh"
#include "3T/P13TOptions.hh"
#include <iosfwd>


// FORWARD REFERENCES

namespace rtt_timestep {
 class ts_manager;
 class field_ts_advisor;
}

namespace XTM {
 class RadiationPhysics;
}

namespace rtt_matprops {
 template<class MT> class MaterialProperties;
}

// DEFINING NAMESPACE

namespace rtt_3T {

 using XTM::RadiationPhysics;

 //===========================================================================//
 // class P13T - 
 // Three Temperature, 3T, P1 and conduction package.
 // 
 //===========================================================================//

 template<class DS>
 class P13T
 {

     // NESTED CLASSES AND TYPEDEFS

   private:

     class P1Coeffs;

   public:

     // Timestep manager typedefs
     
     typedef rtt_timestep::ts_manager ts_manager;
     typedef rtt_timestep::field_ts_advisor field_ts_advisor;

     // Longhand type names.
    
     typedef DS DiffusionSolver;
     typedef typename DiffusionSolver::MeshType MT;
     typedef MT MeshType;

     typedef rtt_matprops::MaterialProperties<MT> MaterialProperties;
    
     // The diffusion solver knows the correct representation for
     // the continuous anb discontinuous face-centered flux fields.

     typedef typename DiffusionSolver::FluxField FluxField;
     typedef typename DiffusionSolver::DiscFluxField DiscFluxField;

     // Miscellaneous shortcut field typedefs from the MT class
    
     typedef typename MeshType::ccsf ccsf;    // cell-centered scalar field

     typedef typename MeshType::ncvsf ncvsf;    // node-centered vector scalar field

     typedef typename MeshType::fcdsf fcdsf;  // face-centered discontinuous s.f.
     typedef typename MeshType::bssf bssf;    // bndry-specified boundary field.

     // The state of the radiation field is passed in and returned
     // in this structure.
    
     struct RadiationStateField {
	 ccsf phi;
	 FluxField F;
	 RadiationStateField(const dsxx::SP<MeshType> &spmesh_)
	     : phi(spmesh_), F(spmesh_)
	 {
	     // empty
	 }
     };
    


     // DATA

   private:
    
     P13TOptions options;            // Specify various solve flags and values
     dsxx::SP<MeshType> spMesh;      // Mesh
     dsxx::SP<ts_manager> spTsManager;           // Timestep Manager
     dsxx::SP<field_ts_advisor> spRadTsAdvisor;  // Radiation Timestep Advisor
     dsxx::SP<field_ts_advisor> spElecTsAdvisor; // Electron Timestep Advisor
     dsxx::SP<field_ts_advisor> spIonTsAdvisor;  // Ion Timestep Advisor
     dsxx::SP<field_ts_advisor> spElecCondTsAdvisor; // Elec conduction advisor
     dsxx::SP<field_ts_advisor> spIonCondTsAdvisor;  // Ion conduction advisor

     // FORBIDDEN METHODS
     
   private:

     // We will not allow copy construction.
     
     P13T(const P13T<DS>& );

     // We will not allow assignment.
     
     P13T& operator=(const P13T& );

   public:

     // CREATORS

     P13T(const P13TOptions &options_, const dsxx::SP<MeshType> &spMesh_);
     P13T(const P13TOptions &options_, const dsxx::SP<MeshType> &spMesh_,
	  dsxx::SP<ts_manager> &spTsManager_);
     ~P13T();

     // MANIPULATORS

     void setOptions(const P13TOptions options_);

     // ACCESSORS

     const dsxx::SP<MeshType> getMesh() const { return spMesh; }
    
     //-----------------------------------------------------------------------//
     // print:
     //     Print itself (for debug mostly)
     //-----------------------------------------------------------------------//

     std::ostream &print(std::ostream &os) const;

     //-----------------------------------------------------------------------//
     // initializeRadiationState:
     //     Initialize the radiation field to Planckian
     //     based on material electron temperatures.
     //-----------------------------------------------------------------------//
    
     void initializeRadiationState(const MaterialProperties &matprops,
				   RadiationStateField &resultsStateField) const;

     //-----------------------------------------------------------------------//
     // solveElectConduction:
     //     Solve for the energy deposition and new temperature due to  
     //     the conduction equation split.
     //-----------------------------------------------------------------------//
    
     void solveElectConduction(ccsf &electronEnergyDeposition,
			       ccsf &Tnp1Electron,
			       DiffusionSolver &solver,
			       double dt,
			       const MaterialProperties &matprops,
			       const bssf &alpha,
			       const bssf &beta,
			       const bssf &bSrc) const;

     //-----------------------------------------------------------------------//
     // solveIonConduction:
     //     Solve for the energy deposition and new temperature due to  
     //     the conduction equation split.
     //-----------------------------------------------------------------------//
    
     void solveIonConduction(ccsf &ionEnergyDeposition,
			     ccsf &Tnp1Ion,
			     DiffusionSolver &solver,
			     double dt,
			     const MaterialProperties &matprops,
			     const bssf &alpha,
			     const bssf &beta,
			     const bssf &bSrc) const;

     //-----------------------------------------------------------------------//
     // solve3T:
     //     Solve for the new radiation field, the electron/ion energy
     //     depositions, and the momentom deposition.
     //-----------------------------------------------------------------------//
    
     void solve3T(RadiationStateField &resultsStateField,
		  ccsf &QEEM,
		  ccsf &REEM,
		  ccsf &electronEnergyDeposition,
		  ccsf &ionEnergyDeposition,
		  ncvsf &momentumDeposition,
		  ccsf &Tnp1Electron,
		  ccsf &Tnp1Ion,
		  DiffusionSolver &solver,
		  double dt,
		  const MaterialProperties &matprops,
		  const ncvsf &velocity,
		  const RadiationStateField &prevStateField,
		  const ccsf &QRad,
		  const ccsf &QElectron,
		  const ccsf &QIon,
		  const bssf &alpha,
		  const bssf &beta,
		  const bssf &bSrc) const;

     // IMPLEMENTATION

   private:
 
     //-----------------------------------------------------------------------//
     // getBhat:
     //    get the 4pi*planckian
     //-----------------------------------------------------------------------//

     void getBhat(ccsf &Bhat, const RadiationPhysics &radPhys,
		  const ccsf &TElectron) const;
    
     //-----------------------------------------------------------------------//
     // getdBhatdT:
     //    get the 4pi*dPlanckiandT
     //-----------------------------------------------------------------------//

     void getdBhatdT(ccsf &dBhatdT, const RadiationPhysics &radPhys,
		     const ccsf &TElectron) const;

     //-----------------------------------------------------------------------//
     // calcNewRadState:
     //     calculate the new radiation state using the previous state,
     //     material properties, and sources.
     //     This solves the coupled radiation, electron, and ion equations
     //     ***without*** the conduction equations.
     //-----------------------------------------------------------------------//
    
     void calcNewRadState(RadiationStateField &resultsStateField,
			  DiffusionSolver &solver,
			  const P1Coeffs &p1coeffs,
			  const bssf &alpha,
			  const bssf &beta,
			  const bssf &bSrc) const;
    
     //-----------------------------------------------------------------------//
     // calcDeltaTElectron:
     //    Calculate the difference between T electron from timestep
     //    n+1 to timestep n+1/2
     //-----------------------------------------------------------------------//

     void calcDeltaTElectron(ccsf &deltaTelectron,
			     double dt,
			     int numGroups, 
			     const P1Coeffs &p1coeffs,
			     const RadiationStateField &radStateField) const;

     //-----------------------------------------------------------------------//
     // calcDeltaTIon:
     //    Calculate the difference between T ion from timestep
     //    n+1 to timestep n+1/2
     //-----------------------------------------------------------------------//

     void calcDeltaTIon(ccsf &deltaTIon,
			double dt,
			const P1Coeffs &p1coeffs,
			const ccsf &deltaTelectron) const;
 };

} // namespace rtt_3T

#include "P1Coeffs.hh"

#endif                          // __3T_P13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13T.hh
//---------------------------------------------------------------------------//
