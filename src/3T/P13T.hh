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

template<class MT, class MaterialProperties, class DiffusionSolver>
class P13T {

    // DATA

  private:
    
    SP<MaterialProperties> spProp;     // Material Props
    P13TOptions options;               // Specify various solve flags and values
    SP<DiffusionSolver> spDiffSolver;  // Which diffusion solver to use


    
    // NESTED CLASSES AND TYPEDEFS

  public:

    // The diffusion solver knows the correct representation for
    // the face-centered current field.

    typedef typename DiffusionSolver::CurrentField CurrentField;

    // The MaterialProperties knows the correct representation for the
    // material state field

    typedef typename MaterialProperties::MaterialStateField MaterialStateField;

    // Miscellaneous shortcut field typedefs from the MT class
    
    typedef typename MT::ccsf ccsf;       // cell-centered scalar field
    typedef typename MT::ncvf ncvf;       // node-centered vector field
    typedef typename MT::fcdcsf fcdsf;    // face-centered discontinuous s.f.

    // The state of the radiation field is passed in and returned
    // in this structure.
    
    struct RadiationStateField {
	ccsf phi;
	CurrentField F;
    };
    


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

    // Initialize the radiation field to Planckian
    // based on material electron temperatures.
    
    void initializeRadiationState(const MaterialStateField &matState,
				  RadiationStateField &resultsStateField) const;

    // Solve for the new radiation field, the electron/ion energy depositions,
    // and the momentom deposition.
    //
    // The P13TOptions object (P13T state variable "options")
    // determines whether this solve is with or without the
    // electron/ion conduction equations.
    
    void solve(double dt,
	       const ccsf QRad,
	       const ccsf QElectron,
	       const ccsf QIon,
	       // **** Boundary Source Field goes here ****
	       const MaterialStateField &matState;
	       const RadiationStateField &prevStateField,
	       RadiationStateField &resultsStateField,
	       ccsf &electronEnergyDeposition,
	       ccsf &ionEnergyDeposition,
	       ncvf &momentumDeposition) const;
};

#endif                          // __3T_P13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P13T.hh
//---------------------------------------------------------------------------//
