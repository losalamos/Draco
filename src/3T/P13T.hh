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

    // The MaterialProperties knows the correct representation for the
    // material state field

    typedef typename MaterialProperties::MaterialStateField MaterialStateField;

    // Miscellaneous shortcut field typedefs from the MT class
    
    typedef typename MeshType::ccsf ccsf;    // cell-centered scalar field
    typedef typename MeshType::ncvf ncvf;    // node-centered vector field
    typedef typename MeshType::fcdcsf fcdsf; // face-centered discontinuous s.f.
    typedef typename MeshType::bsbf bsbf;    // bndry-specified boundary field.

    // The state of the radiation field is passed in and returned
    // in this structure.
    
    struct RadiationStateField {
	ccsf phi;
	FluxField F;
    };
    


    // DATA

  private:
    
    P13TOptions options;               // Specify various solve flags and values
    SP<MaterialProperties> spProp;     // Material Props
    SP<DiffusionSolver> spDiffSolver;  // Which diffusion solver to use


    
  public:

    // CREATORS

    P13T(const P13TOptions &options_,
	 const SP<MaterialProperties> &spProp_,
	 const SP<DiffusionSolver> &spDiffSolver_);
    P13T(const P13T& );
    ~P13T();

    // MANIPULATORS

    P13T& operator=(const P13T& );
    void setMaterialProperties(const SP<MaterialProperties> &spProp_);
    void setOptions(const P13TOptions options_);
    void setDiffSolver(const SP<DiffusionSolver> &spDiffSolver_);

    // ACCESSORS

    //------------------------------------------------------------------------//
    // initializeRadiationState:
    //     Initialize the radiation field to Planckian
    //     based on material electron temperatures.
    //------------------------------------------------------------------------//
    
    void initializeRadiationState(const MaterialStateField &matState,
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
	       const MaterialStateField &matState,
	       const RadiationStateField &prevStateField,
	       const ccsf QRad,
	       const ccsf QElectron,
	       const ccsf QIon,
	       const bsbf boundary,
	       RadiationStateField &resultsStateField,
	       ccsf &electronEnergyDeposition,
	       ccsf &ionEnergyDeposition,
	       ncvf &momentumDeposition) const;

    // IMPLEMENTATION

  private:
    
    //------------------------------------------------------------------------//
    // clacNewRadState:
    //     calculate the new radiation state using the previous state,
    //     material properties, and sources.
    //     This solves the coupled radiation, electron, and ion equations
    //     ***without*** the conduction equations.
    //------------------------------------------------------------------------//
    
    void calcNewRadState(double dt,
			 int groupNo,
			 const MaterialStateField &matState,
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
		      const MaterialStateField &matState,
		      const RadiationStateField &prevStateField,
		      const ccsf QRad,
		      const ccsf QElectron,
		      const ccsf QIon,
		      const ccsf TElectron,
		      const ccsf TIon,
		      fcdsf &D,
		      FluxField &Fprime,
		      ccsf &sigmaAbsBar,
		      ccsf &QRadBar) const;
};

#endif                          // __3T_P13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13T.hh
//---------------------------------------------------------------------------//
