//----------------------------------*-C++-*----------------------------------//
// P1Coeffs.hh
// Randy M. Roberts
// Thu Nov  5 13:14:06 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_P1Coeffs_hh__
#define __3T_P1Coeffs_hh__

#include "3T/P13T.hh"

#include "ds++/SP.hh"

namespace rtt_3T
{
 //===========================================================================//
 // class P1Coeffs - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class DS>
 class P13T<DS>::P1Coeffs
 {

     // NESTED CLASSES AND TYPEDEFS

     typedef typename P13T<DS>::fcdsf fcdsf;
     typedef typename P13T<DS>::ccsf ccsf;
     typedef typename P13T<DS>::DiscFluxField DiscFluxField;

     // DATA
    
     dsxx::SP<fcdsf> spD;
     dsxx::SP<DiscFluxField> spFprime;
     dsxx::SP<ccsf> spSigmaAbsBar;
     dsxx::SP<ccsf> spQEEM;
     dsxx::SP<ccsf> spREEM;
     dsxx::SP<ccsf> spQRadBar;
     dsxx::SP<ccsf> spQElecStar;
     dsxx::SP<ccsf> spCvStar;
     dsxx::SP<ccsf> spNu;

     const P13T<DS> &p13T;
     const dsxx::SP<MT> spMesh;
     double dt;
     int    groupNo;
     const P13TOptions &options;
     const MaterialProperties &matprops;
     const ncvsf &velocity;
     const RadiationStateField &prevStateField;
     const ccsf &QRad;
     const ccsf &QElectron;
     const ccsf &QIon;
     const ccsf &TElectron;
     const ccsf &TIon;

   public:

     // CREATORS
    
     P1Coeffs(const P13T<DS> &p13T_,
	      double dt_,
	      int groupNo_,
	      const P13TOptions &options_,
	      const MaterialProperties &matprops_,
	      const ncvsf &velocity_,
	      const RadiationStateField &prevStateField_,
	      const ccsf &QRad_,
	      const ccsf &QElectron_,
	      const ccsf &QIon_,
	      const ccsf &TElectron_,
	      const ccsf &TIon_);

     // MANIPULATORS
    
     // ACCESSORS

     const fcdsf &D() const { return *spD; }
     const DiscFluxField &Fprime() const { return *spFprime; }
     const ccsf &sigmaAbsBar() const { return *spSigmaAbsBar; }
     const ccsf &QEEM() const { return *spQEEM; }
     const ccsf &REEM() const { return *spREEM; }
     const ccsf &QRadBar() const { return *spQRadBar; }
     const ccsf &QElecStar() const { return *spQElecStar; }
     const ccsf &CvStar() const { return *spCvStar; }
     const ccsf &nu() const { return *spNu; }

     double get_dt() const { return dt; }
     const MaterialProperties &getMatprops() const { return matprops; }
     const ccsf &getQRad() const { return QRad; }
     const ccsf &getQElectron() const { return QElectron; }
     const ccsf &getQIon() const { return QIon; }
     const ccsf &getTElectron() const { return TElectron; }
     const ccsf &getTIon() const { return TIon; }

   private:
    
     // IMPLEMENTATION
     
     //-----------------------------------------------------------------------//
     // calcP1Coeffs:
     //     Calculate the coefficients, e.g. diffusion and removal, and
     //     source terms for solving the P1 equation.
     //-----------------------------------------------------------------------//

     void calcP1Coeffs();

     //-----------------------------------------------------------------------//
     // calcStarredFieldsAndNu:
     //    Calculate Qe*, Cv*, and nu.
     //    These are needed to calculate other coefficients
     //    and delta temperatures.
     //-----------------------------------------------------------------------//

     void calcStarredFieldsAndNu();
    
     //-----------------------------------------------------------------------//
     // calcStarredFields:
     //    Calculate Qe*, Cv*, but not nu.
     //    These are needed to calculate other coefficients
     //    and delta temperatures.
     //-----------------------------------------------------------------------//

     void calcStarredFields();    

 };

} // end namespace rtt_3T

#endif                          // __3T_P1Coeffs_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/P1Coeffs.hh
//---------------------------------------------------------------------------//
