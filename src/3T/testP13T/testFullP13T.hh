//----------------------------------*-C++-*----------------------------------//
// testFullP13T.hh
// Randy M. Roberts
// Sat May 30 15:09:58 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_testFullP13T_hh__
#define __3T_testP13T_testFullP13T_hh__

#include "mesh/Mesh_XYZ.hh"
#include "testFullP13T_DB.hh"
#include "ds++/SP.hh"
#include "3T/P13T.hh"
// #include "3T/Diffusion_P1.hh"
#include "P1Diffusion/P1Diffusion.hh"
#include "P1Diffusion/SolverP1Diff.hh"

#define MARSHAK_MATPROPS

#include "matprops/MarshakMaterialProps.hh"
#include "matprops/InterpedMaterialProps.hh"

#include <string>

// FORWARD REFERENCES

namespace rtt_timestep {
 class ts_manager;
}

namespace rtt_matprops {
 template<class MT> class TempMapper;
}

namespace XTM {
 
 //===========================================================================//
 // class testFullP13T - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 class testFullP13T
 {

     // NESTED CLASSES AND TYPEDEFS

   public:
    
     typedef Mesh_XYZ MT;

     typedef rtt_matprops::MarshakMaterialProps MarshakMaterialProps;
     typedef rtt_matprops::InterpedMaterialProps InterpedMaterialProps;

#ifdef MARSHAK_MATPROPS
     typedef MarshakMaterialProps MP;
#else
     typedef InterpedMaterialProps MP;
#endif
     
     // typedef Diffusion_P1<MT> DS;

     typedef rtt_P1Diffusion::SolverP1Diff<MT> MS;
     typedef MS MatrixSolver;
     typedef rtt_P1Diffusion::P1Diffusion<MT, MS> DS;

     typedef MT::ccsf ccsf;
     typedef MT::ccif ccif;
     typedef MT::fcdsf fcdsf;
     typedef MT::fcdif fcdif;
     typedef MT::bssf bssf;
     
     typedef P13T<MT,MP,DS> P13T;

     typedef P13T::RadiationStateField RadiationStateField;
     
     typedef MP::MaterialStateField<MT::ccsf> MatStateCC;
     typedef MP::MaterialStateField<MT::fcdsf> MatStateFC;
     
     // DATA

   private:
    
     Units units;
     dsxx::SP<MP> spMatProp;
     dsxx::SP<rtt_matprops::TempMapper<MT> > spTempMapper;

     dsxx::SP<MT> spMesh;
     dsxx::SP<P13T> spP13T;

     dsxx::SP<rtt_timestep::ts_manager> spTsManager;

     mutable dsxx::SP<DS> spDiffSolver;
     
     testFullP13T_DB pdb;
     Diffusion_DB diffdb;
     pcg_DB pcg_db;
    
   public:

     // CREATORS
    
     testFullP13T(const std::string &infile);
     ~testFullP13T();

     // MANIPULATORS
    
     // ACCESSORS

     void run() const;

   private:
    
     // IMPLEMENTATION

     void timestep(double &time, double &dt, int &cycle,
		   MatStateCC &matStateCC, MatStateFC &matStateFC,
		   RadiationStateField &radState,
		   ccsf &electEnergyDep, ccsf &ionEnergyDep,
		   const ccsf &QRad, const ccsf &QElectron, const ccsf &QIon,
		   const bssf &alpha, const bssf &beta, const bssf &bSrc) const;

     void getMatProp();
     void getMatProp(dsxx::SP<MarshakMaterialProps> &spMatProp_) const;
     void getMatProp(dsxx::SP<InterpedMaterialProps> &spMatProp_) const;
    
     MatStateCC getMatStateCC(const ccsf &TElec, const ccsf &TIon,
			      const ccsf &density, const ccif &matid) const;

     MatStateFC getMatStateFC(const ccsf &TElec, const ccsf &TIon,
			      const ccsf &density, const ccif &matid) const;

     void gmvDump(const RadiationStateField &radState, const ccsf &TElec,
		  const ccsf &TIon, int cycle, double time) const;

     void setBoundary(bssf &alpha, bssf &beta, bssf &bSrc) const;
     
     void postProcess(const RadiationStateField &radState,
		      const RadiationStateField &newRadState,
		      const MatStateCC &matStateCC,
		      const MatStateCC &newMatStateCC,
		      const ccsf &electEnergyDep,
		      const ccsf &ionEnergyDep,
		      const ccsf &QRad,
		      const ccsf &QElectron,
		      const ccsf &QIon,
		      const ccsf &QEEM,
		      const ccsf &REEM,
		      double dt) const;
     
     double calcLeakage(const RadiationStateField &radstate) const;
     
 };

}  // namespace XTM

#endif                          // __3T_testP13T_testFullP13T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/testFullP13T.hh
//---------------------------------------------------------------------------//
