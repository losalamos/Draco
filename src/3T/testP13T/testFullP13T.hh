//----------------------------------*-C++-*----------------------------------//
// testFullP13T.hh
// Randy M. Roberts
// Sat May 30 15:09:58 1998
//---------------------------------------------------------------------------//
// @> Test facility for the P1 3T radiation solver.
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
#include "matprops/MultiMatCellMatProps.hh"

#include <string>

// FORWARD REFERENCES

namespace rtt_timestep {
    class ts_manager;
    class fixed_ts_advisor;
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
 
 template<class UMCMP>
 class testFullP13T
 {
     
     // NESTED CLASSES AND TYPEDEFS
     
   public:
    
     typedef Mesh_XYZ MT;

     typedef rtt_matprops::MarshakMaterialProps  MarshakMaterialProps;
     typedef rtt_matprops::InterpedMaterialProps InterpedMaterialProps;
     typedef rtt_matprops::MultiMatCellMatProps<UMCMP> MP;
     
     typedef MarshakMaterialProps  MP_MRSHK;
     typedef InterpedMaterialProps MP_INTRP;

     // typedef Diffusion_P1<MT> DS;

     typedef rtt_P1Diffusion::SolverP1Diff<MT> MS;
     typedef MS MatrixSolver;
     typedef rtt_P1Diffusion::P1Diffusion<MT, MS> DS;

     typedef MT::ccsf  ccsf;
     typedef MT::ccif  ccif;
     typedef MT::fcdsf fcdsf;
     typedef MT::fcdif fcdif;
     typedef MT::bssf  bssf;

     typedef MT::cctf  < std::vector<double  > > ccvsf;
     typedef MT::fcdtf < std::vector<double  > > fcdvsf;
     typedef MT::cctf  < std::vector<int     > > ccvif;
     typedef MT::fcdtf < std::vector<int     > > fcdvif;
     
     typedef typename MP:: template MaterialStateField<ccsf, ccvsf, ccvif> MatStateCC;
     typedef typename MP:: template MaterialStateField<fcdsf, fcdvsf, fcdvif> MatStateFC;

     typedef P13T<MT,MatStateCC,MatStateFC,DS> P13T;
     typedef typename P13T::RadiationStateField RadiationStateField;

     // DATA

   private:
    
     Units units;

     dsxx::SP<MP>       spMatProp;
     dsxx::SP<UMCMP>    spUMatProp;
     dsxx::SP<rtt_matprops::TempMapper<MT> > spTempMapper;

     dsxx::SP<MT> spMesh;
     dsxx::SP<P13T> spP13T;

     dsxx::SP<rtt_timestep::ts_manager> spTsManager;
     dsxx::SP<rtt_timestep::fixed_ts_advisor> spTsCurrent;

     mutable dsxx::SP<DS> spDiffSolver;
     
     testFullP13T_DB tdb;
     Diffusion_DB diffdb;
     pcg_DB pcg_db;
    
   public:

     // CREATORS
    
     testFullP13T(const testFullP13T_DB &tdb_,
		  const Diffusion_DB &diffdb,
		  const Mesh_DB &mdb,
		  const pcg_DB &pgc_db_);

     ~testFullP13T();

     // MANIPULATORS
    
     // ACCESSORS

     void run() const;
     const Units &getUnits() const { return units; }

   private:
    
     // IMPLEMENTATION


     void timestep(double &time, double &dt, int &cycle,
		   MatStateCC &matStateCC, MatStateFC &matStateFC,
		   RadiationStateField &radState,
		   ccsf &electEnergyDep, ccsf &ionEnergyDep,
		   const ccsf &QRad, const ccsf &QElectron, const ccsf &QIon,
		   const bssf &alpha, const bssf &beta, const bssf &bSrc) const;

     void getMatProp();
    
     MatStateCC getMatStateCC(const ccvsf &TElec, const ccvsf &TIon,
			      const ccvsf &density, const ccvsf &VolFrac,
			      const ccvif &matid) const;

     MatStateFC getMatStateFC(const MatStateCC &msfcc) const;

     void gmvDump(const RadiationStateField &radState, const ccsf &TElec,
		  const ccsf &TIon, int dumpno, int cycle, double time) const;

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
