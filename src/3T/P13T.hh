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

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

// FORWARD REFERENCES

class RadiationPhysics;

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
    typedef typename MeshType::ncvf ncvf;    // node-centered vector field
    typedef typename MeshType::fcdsf fcdsf;  // face-centered discontinuous s.f.
    typedef typename MeshType::bsbf bsbf;    // bndry-specified boundary field.

    // The MaterialProperties knows the correct representation for the
    // material state field

    typedef typename MaterialProperties::MaterialStateField<ccsf>
            CCMaterialStateField;

    typedef typename MaterialProperties::MaterialStateField<fcdsf>
            FCMaterialStateField;

    // The state of the radiation field is passed in and returned
    // in this structure.
    
    struct RadiationStateField {
	ccsf phi;
	FluxField F;
	RadiationStateField(const SP<MeshType> &spmesh_)
	    : phi(spmesh_), F(spmesh_)
	{
	    // empty
	}
    };
    


    // DATA

  private:
    
    P13TOptions options;               // Specify various solve flags and values
    SP<MeshType> spMesh;               // Mesh
    SP<DiffusionSolver> spDiffSolver;  // Which diffusion solver to use


    
  public:

    // CREATORS

    P13T(const P13TOptions &options_,
	 const SP<DS> &spDiffSolver_);
    P13T(const P13T<MT,MP,DS>& );
    ~P13T();

    // MANIPULATORS

    P13T& operator=(const P13T& );
    void setOptions(const P13TOptions options_);
    void setDiffSolver(const SP<DiffusionSolver> &spDiffSolver_);

    // ACCESSORS

    const SP<MeshType> getMesh() const { return spMesh; }
    
    //------------------------------------------------------------------------//
    // print:
    //     Print itself (for debug mostly)
    //------------------------------------------------------------------------//

    std::ostream &print(std::ostream &os) const;

    //------------------------------------------------------------------------//
    // initializeRadiationState:
    //     Initialize the radiation field to Planckian
    //     based on material electron temperatures.
    //------------------------------------------------------------------------//
    
    void initializeRadiationState(const CCMaterialStateField &matStateCC,
				  RadiationStateField &resultsStateField) const;

    //------------------------------------------------------------------------//
    // solve:
    //     Solve for the new radiation field, the electron/ion energy
    //     depositions, and the momentom deposition.
    //
    //     The P13TOptions object (P13T state variable "options")
    //     determines whether this solve is with or without the
    //     electron/ion conduction equations.
    //------------------------------------------------------------------------//
    
    void solve(double dt,
	       const CCMaterialStateField &matStateCC,
	       const FCMaterialStateField &matStateFC,
	       const RadiationStateField &prevStateField,
	       const ccsf QRad,
	       const ccsf QElectron,
	       const ccsf QIon,
	       const bsbf boundary,
	       RadiationStateField &resultsStateField,
	       ccsf &electronEnergyDeposition,
	       ccsf &ionEnergyDeposition,
	       ncvf &momentumDeposition,
	       ccsf &Tnp1Electron,
	       ccsf &Tnp1Ion) const;

    // IMPLEMENTATION

  private:
 
    //------------------------------------------------------------------------//
    // getBhat:
    //    get the 4pi*planckian
    //------------------------------------------------------------------------//

    void getBhat(const RadiationPhysics &radPhys,
		 const ccsf &TElectron, ccsf &Bhat) const;
    
    //------------------------------------------------------------------------//
    // getdBhatdT:
    //    get the 4pi*dPlanckiandT
    //------------------------------------------------------------------------//

    void getdBhatdT(const RadiationPhysics &radPhys,
		    const ccsf &TElectron, ccsf &dBhatdT) const;

    //------------------------------------------------------------------------//
    // calcNewRadState:
    //     calculate the new radiation state using the previous state,
    //     material properties, and sources.
    //     This solves the coupled radiation, electron, and ion equations
    //     ***without*** the conduction equations.
    //------------------------------------------------------------------------//
    
    void calcNewRadState(double dt,
			 int groupNo,
			 const CCMaterialStateField &matStateCC,
			 const FCMaterialStateField &matStateFC,
			 const RadiationStateField &prevStateField,
			 const ccsf QRad,
			 const ccsf QElectron,
			 const ccsf QIon,
			 const ccsf TElectron,
			 const ccsf TIon,
			 const bsbf boundary,
			 RadiationStateField &resultsStateField) const;
    
    //------------------------------------------------------------------------//
    // calcP1Coeffs:
    //     Calculate the coefficients, e.g. diffusion and removal, and
    //     source terms for solving the P1 equation.
    //------------------------------------------------------------------------//

    void calcP1Coeffs(double dt,
		      int groupNo,
		      const CCMaterialStateField &matStateCC,
		      const FCMaterialStateField &matStateFC,
		      const RadiationStateField &prevStateField,
		      const ccsf &QRad,
		      const ccsf &QElectron,
		      const ccsf &QIon,
		      const ccsf &TElectron,
		      const ccsf &TIon,
		      fcdsf &D,
		      DiscFluxField &Fprime,
		      ccsf &sigmaAbsBar,
		      ccsf &QRadBar) const;

    //------------------------------------------------------------------------//
    // calcStarredFields:
    //    Calculate Qe*, Cv*, and nu.
    //    These are needed to calculate other coefficients
    //    and delta temperatures.
    //------------------------------------------------------------------------//

    void calcStarredFields(double dt,
			   int groupNo,
			   const CCMaterialStateField &matStateCC,
			   const RadiationPhysics &radPhys,
			   const ccsf &QElectron,
			   const ccsf &QIon,
			   const ccsf &TElectron,
			   const ccsf &TIon,
			   const ccsf &sigmaEmission,
			   ccsf &QElecStar,
			   ccsf &CvStar,
			   ccsf &nu) const;
    
    //------------------------------------------------------------------------//
    // calcStarredFields:
    //    Calculate Qe*, Cv*, but not nu.
    //    These are needed to calculate other coefficients
    //    and delta temperatures.
    //------------------------------------------------------------------------//

    void calcStarredFields(double dt,
			   int groupNo,
			   const CCMaterialStateField &matStateCC,
			   const RadiationPhysics &radPhys,
			   const ccsf &QElectron,
			   const ccsf &QIon,
			   const ccsf &TElectron,
			   const ccsf &TIon,
			   const ccsf &sigmaEmission,
			   ccsf &QElecStar,
			   ccsf &CvStar) const;
    
    //------------------------------------------------------------------------//
    // calcDeltaTElectron:
    //    Calculate the difference between T electron from timestep
    //    n+1 to timestep n+1/2
    //------------------------------------------------------------------------//

    void calcDeltaTElectron(double dt,
			    int numGroups, 
			    const CCMaterialStateField &matStateCC, 
			    const RadiationStateField &prevStateField, 
			    const ccsf &QElectron, 
			    const ccsf &QIon,
			    const ccsf &TElectron,
			    const ccsf &TIon,
			    const RadiationStateField &resultsStateField, 
			    ccsf &deltaTelectron) const;

    //------------------------------------------------------------------------//
    // calcDeltaTIon:
    //    Calculate the difference between T ion from timestep
    //    n+1 to timestep n+1/2
    //------------------------------------------------------------------------//

    void calcDeltaTIon(double dt,
		       const CCMaterialStateField &matStateCC, 
		       const RadiationStateField &prevStateField, 
		       const ccsf &QIon,
		       const ccsf &TElectron,
		       const ccsf &TIon,
		       const ccsf &deltaTelectron,
		       ccsf &deltaTIon) const;
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

END_NS_XTM  // namespace XTM

#endif                          // __3T_P13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13T.hh
//---------------------------------------------------------------------------//
