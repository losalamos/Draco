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

// DEFINING NAMESPACE

namespace XTM {

 //===========================================================================//
 // class P13T - 
 // Three Temperature, 3T, P1 and conduction package.
 // 
 //===========================================================================//

 template<class MT, class MP, class DS>
 class P13T
 {

     // NESTED CLASSES AND TYPEDEFS

   public:

     // Timestep manager typedefs
     
     typedef rtt_timestep::ts_manager ts_manager;
     typedef rtt_timestep::field_ts_advisor field_ts_advisor;

     // Longhand type names.
    
     typedef MT MeshType;
     typedef DS DiffusionSolver;
     typedef MP MaterialProperties;
    
     // The diffusion solver knows the correct representation for
     // the continuous anb discontinuous face-centered flux fields.

     typedef typename DiffusionSolver::FluxField FluxField;
     typedef typename DiffusionSolver::DiscFluxField DiscFluxField;

     // Miscellaneous shortcut field typedefs from the MT class
    
     typedef typename MeshType::ccsf ccsf;    // cell-centered scalar field
#if 0
     typedef typename MeshType::ncvf ncvf;    // node-centered vector field
#endif
     typedef typename MeshType::fcdsf fcdsf;  // face-centered discontinuous s.f.
     typedef typename MeshType::bssf bssf;    // bndry-specified boundary field.

     // The MaterialProperties knows the correct representation for the
     // material state field

     typedef typename MaterialProperties::template MaterialStateField<ccsf>
     CCMaterialStateField;

     typedef typename MaterialProperties::template MaterialStateField<fcdsf>
     FCMaterialStateField;

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
     
     P13T(const P13T<MT,MP,DS>& );

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
    
     void initializeRadiationState(const CCMaterialStateField &matStateCC,
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
			       const CCMaterialStateField &matStateCC,
			       const FCMaterialStateField &matStateFC) const;

     //-----------------------------------------------------------------------//
     // solveIonConduction:
     //     Solve for the energy deposition and new temperature due to  
     //     the conduction equation split.
     //-----------------------------------------------------------------------//
    
     void solveIonConduction(ccsf &ionEnergyDeposition,
			     ccsf &Tnp1Ion,
			     DiffusionSolver &solver,
			     double dt,
			     const CCMaterialStateField &matStateCC,
			     const FCMaterialStateField &matStateFC) const;

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
#if 0
		  ncvf &momentumDeposition,
#endif
		  ccsf &Tnp1Electron,
		  ccsf &Tnp1Ion,
		  DiffusionSolver &solver,
		  double dt,
		  const CCMaterialStateField &matStateCC,
		  const FCMaterialStateField &matStateFC,
		  const RadiationStateField &prevStateField,
		  const ccsf &QRad,
		  const ccsf &QElectron,
		  const ccsf &QIon,
		  const bssf &boundary) const;

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
			  ccsf &QEEM,
			  ccsf &REEM,
			  DiffusionSolver &solver,
			  double dt,
			  int groupNo,
			  const CCMaterialStateField &matStateCC,
			  const FCMaterialStateField &matStateFC,
			  const RadiationStateField &prevStateField,
			  const ccsf &QRad,
			  const ccsf &QElectron,
			  const ccsf &QIon,
			  const ccsf &TElectron,
			  const ccsf &TIon,
			  const bssf &boundary) const;
    
     //-----------------------------------------------------------------------//
     // calcP1Coeffs:
     //     Calculate the coefficients, e.g. diffusion and removal, and
     //     source terms for solving the P1 equation.
     //-----------------------------------------------------------------------//

     void calcP1Coeffs(fcdsf &D,
		       DiscFluxField &Fprime,
		       ccsf &sigmaAbsBar,
		       ccsf &QEEM,
		       ccsf &REEM,
		       ccsf &QRadBar,
		       double dt,
		       int groupNo,
		       const CCMaterialStateField &matStateCC,
		       const FCMaterialStateField &matStateFC,
		       const RadiationStateField &prevStateField,
		       const ccsf &QRad,
		       const ccsf &QElectron,
		       const ccsf &QIon,
		       const ccsf &TElectron,
		       const ccsf &TIon) const;

     //-----------------------------------------------------------------------//
     // calcStarredFields:
     //    Calculate Qe*, Cv*, and nu.
     //    These are needed to calculate other coefficients
     //    and delta temperatures.
     //-----------------------------------------------------------------------//

     void calcStarredFields(ccsf &QElecStar,
			    ccsf &CvStar,
			    ccsf &nu,
			    double dt,
			    int groupNo,
			    const CCMaterialStateField &matStateCC,
			    const RadiationPhysics &radPhys,
			    const ccsf &QElectron,
			    const ccsf &QIon,
			    const ccsf &TElectron,
			    const ccsf &TIon,
			    const ccsf &sigmaEmission) const;
    
     //-----------------------------------------------------------------------//
     // calcStarredFields:
     //    Calculate Qe*, Cv*, but not nu.
     //    These are needed to calculate other coefficients
     //    and delta temperatures.
     //-----------------------------------------------------------------------//

     void calcStarredFields(ccsf &QElecStar,
			    ccsf &CvStar,
			    double dt,
			    int groupNo,
			    const CCMaterialStateField &matStateCC,
			    const RadiationPhysics &radPhys,
			    const ccsf &QElectron,
			    const ccsf &QIon,
			    const ccsf &TElectron,
			    const ccsf &TIon) const;
    
     //-----------------------------------------------------------------------//
     // calcDeltaTElectron:
     //    Calculate the difference between T electron from timestep
     //    n+1 to timestep n+1/2
     //-----------------------------------------------------------------------//

     void calcDeltaTElectron(ccsf &deltaTelectron,
			     double dt,
			     int numGroups, 
			     const CCMaterialStateField &matStateCC, 
			     const RadiationStateField &prevStateField, 
			     const ccsf &QElectron, 
			     const ccsf &QIon,
			     const ccsf &TElectron,
			     const ccsf &TIon,
			     const RadiationStateField &resultsStateField) const;

     //-----------------------------------------------------------------------//
     // calcDeltaTIon:
     //    Calculate the difference between T ion from timestep
     //    n+1 to timestep n+1/2
     //-----------------------------------------------------------------------//

     void calcDeltaTIon(ccsf &deltaTIon,
			double dt,
			const CCMaterialStateField &matStateCC, 
			const RadiationStateField &prevStateField, 
			const ccsf &QIon,
			const ccsf &TElectron,
			const ccsf &TIon,
			const ccsf &deltaTelectron) const;
 };

 //---------------------------------------------------------------------------//
 // operator<<:
 //    A convenience function to print a P13T
 //    (mostly for debug purposes)
 //---------------------------------------------------------------------------//

 template<class MT, class MP, class DS>
 inline std::ostream &operator<<(std::ostream &os, const P13T<MT, MP, DS> &rhs)
{
    return rhs.print(os);
}

} // namespace XTM

#endif                          // __3T_P13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13T.hh
//---------------------------------------------------------------------------//
