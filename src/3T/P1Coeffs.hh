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

     // DATA
    
     const P13T<DS> &p13T;
     const dsxx::SP<MT> spMesh;

     fcdsf xD;
     DiscFluxField xFprime;
     ccsf xSigmaAbsBar;
     ccsf xQEEM;
     ccsf xREEM;
     ccsf xQRadBar;
     ccsf xQElecStar;
     ccsf xCvStar;
     ccsf xNu;
     ccsf xXiTilde;
     ccsf xXiBracket;
     fcdsf xXiOmegaBracket;

     double dt;
     int    groupNo;
     const DiffusionSolver &solver;
     const P13TOptions &options;
     const MaterialProperties &matprops;
#ifdef P13T_MOMENTUM_DEPOSITION
     const ncvsf &velocity;
#endif
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
              const DiffusionSolver &solver_,
	      const MaterialProperties &matprops_,
#ifdef P13T_MOMENTUM_DEPOSITION
	      const ncvsf &velocity_,
#endif	      
	      const RadiationStateField &prevStateField_,
	      const ccsf &QRad_,
	      const ccsf &QElectron_,
	      const ccsf &QIon_,
	      const ccsf &TElectron_,
	      const ccsf &TIon_);

     // MANIPULATORS
    
     // ACCESSORS

     const fcdsf &D() const { return xD; }
     const DiscFluxField &Fprime() const { return xFprime; }
     const ccsf &sigmaAbsBar() const { return xSigmaAbsBar; }
     const ccsf &QEEM() const { return xQEEM; }
     const ccsf &REEM() const { return xREEM; }
     const ccsf &QRadBar() const { return xQRadBar; }
     const ccsf &QElecStar() const { return xQElecStar; }
     const ccsf &CvStar() const { return xCvStar; }
     const ccsf &nu() const { return xNu; }
     const ccsf &xiTilde() const { return xXiTilde; }
     const ccsf &xiBracket() const { return xXiBracket; }
     const fcdsf &xiOmegaBracket() const { return xXiOmegaBracket; }

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
     // calcVelocityCorrections:
     //    Calculate the velocity correction terms.
     //-----------------------------------------------------------------------//

     void calcVelocityCorrections();

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
