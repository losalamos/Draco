//----------------------------------*-C++-*----------------------------------//
// P1Coeffs.t.cc
// Randy M. Roberts
// Thu Nov  5 13:14:06 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/P1Coeffs.hh"

namespace rtt_3T
{

 //---------------------------------------------------------------------------//
 // P1Coeffs:
 //     The constructur for the P1Coeffs object.
 //     We will allocate memory for our smart pointers,
 //     and then call calcP1Coeffs to load them up.
 //---------------------------------------------------------------------------//
 
 template<class DS>
 P13T<DS>::P1Coeffs::P1Coeffs(const P13T<DS> &p13T_,
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
			      const ccsf &TIon_)
     : p13T(p13T_), spMesh(p13T_.getMesh()),
       dt(dt_), groupNo(groupNo_), options(options_),
       matprops(matprops_), velocity(velocity_),
       prevStateField(prevStateField_),
       QRad(QRad_), QElectron(QElectron_), QIon(QIon_),
       TElectron(TElectron_), TIon(TIon_)
 {
     spQEEM = new ccsf(spMesh);
     spREEM = new ccsf(spMesh);
     spQRadBar = new ccsf(spMesh);
     spQElecStar = new ccsf(spMesh);
     spCvStar = new ccsf(spMesh);
     spNu = new ccsf(spMesh);
     spSigmaAbsBar = new ccsf(spMesh);
     spD = new fcdsf(spMesh);
     spFprime = new DiscFluxField(spMesh);
     
     calcP1Coeffs();
 }

 //---------------------------------------------------------------------------//
 // calcP1Coeffs:
 //     Calculate the coefficients, e.g. diffusion and removal, and
 //     source terms for solving the P1 equation.
 //---------------------------------------------------------------------------//

 template<class DS>
 void P13T<DS>::P1Coeffs::calcP1Coeffs()
 {
     // Set the radiation physics to the given units.
    
     const RadiationPhysics radPhys(matprops.getUnits());

     // set up some needed scalars, like tau

     double c = radPhys.getLightSpeed();
     double tau = 1.0/(c*dt);

    // If this is diffusion, then tauP1 = 0.0.
    // If this is full P1 then tauP1 = tau.
    
     double tauP1 = tau * options.getP1TauMultiplier();

     // Ask the material properties for sigma total.
     // It is the material properties' responsibility to do
     // any averaging of temperatures, etc. to achieve the correct
     // resulting sigmaTotal.

     fcdsf sigmaTotal(spMesh);
     matprops.getSigmaTotal(groupNo, sigmaTotal);

     //
     // We can now calculate the results
     //

     ccsf sigmaAbs(spMesh);
     matprops.getSigmaAbsorption(groupNo, sigmaAbs);
    
     // Allocate and calculate the diffusion constant.

     Assert(spD);
     
     *spD = (1.0/3.0) / (sigmaTotal + tauP1);

     // Make shorthand references.

     Assert(spQEEM);
     Assert(spREEM);
     Assert(spQRadBar);
     Assert(spQElecStar);
     Assert(spCvStar);
     Assert(spNu);
     
     ccsf &QEEM = *spQEEM;
     ccsf &REEM = *spREEM;
     ccsf &QRadBar = *spQRadBar;
     const ccsf &QElecStar = *spQElecStar;
     const ccsf &CvStar = *spCvStar;
     const ccsf &nu = *spNu;

     if (options.getIsCoupledMaterial())
     {
	 calcStarredFieldsAndNu();

	 cerr << "QElecStar: " << *QElecStar.begin() << endl;
	 cerr << "CvStar: " << *CvStar.begin() << endl;
	 cerr << "nu: " << *nu.begin() << endl;
    
	 // Calculate modified sigma absorption

	 Assert(spSigmaAbsBar);
	 
	 *spSigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

	 // Calculate the emmissive removal coefficient
	
	 REEM = -1.0 * nu * sigmaAbs;
	
	 // We need the Bhat.

	 ccsf Bhat(spMesh);
	 p13T.getBhat(Bhat, radPhys, TElectron);

	 // Get sigma emission
    
	 ccsf sigmaEmission(spMesh);
	 matprops.getSigmaEmission(groupNo, sigmaEmission);

	 // Calculate the emmissive source term.
	
	 QEEM = (1.0 - nu)*sigmaEmission*Bhat + nu*QElecStar;

	 // Calculated modified radiation source

	 QRadBar = tau*prevStateField.phi + QEEM + QRad;
     }
     else
     {
	 Assert(spSigmaAbsBar);
	 *spSigmaAbsBar = sigmaAbs + tau;
	 REEM = 0.0;
	 QEEM = 0.0;
	 QRadBar = tau*prevStateField.phi + QRad;
     }

     // Calculate the "telegraph" term to the P1 equation.

     Assert(spFprime);
     *spFprime = tauP1*prevStateField.F / (sigmaTotal + tauP1);

 }

 //------------------------------------------------------------------------//
 // calcStarredFieldsAndNu:
 //    Calculate Qe*, Cv*, and nu.
 //    These are needed to calculate other coefficients
 //    and delta temperatures.
 //------------------------------------------------------------------------//

 template<class DS>
 void P13T<DS>::P1Coeffs::calcStarredFieldsAndNu()
 {
     // Set the radiation physics to the given units.
    
     const RadiationPhysics radPhys(matprops.getUnits());

     // Calculate Qe* and Cv*.
     // We will then calculate nu ourself.
    
     calcStarredFields();

    // Calculate the 4pi*Planckian's temperature derivative.

     ccsf dBhatdT(spMesh);
     p13T.getdBhatdT(dBhatdT, radPhys, TElectron);
    
    // Get sigma emission
    
     ccsf sigmaEmission(spMesh);
     matprops.getSigmaEmission(groupNo, sigmaEmission);

    // Make shorthand references.

     Assert(spCvStar);
     Assert(spNu);
     
     const ccsf &CvStar = *spCvStar;
     ccsf &nu = *spNu;

    // Calculate the "nu" used in the 3T modification of sigmaAbs
    
     nu = dt * sigmaEmission * dBhatdT /
	 (CvStar + dt * sigmaEmission * dBhatdT);
 }

 //------------------------------------------------------------------------//
 // calcStarredFields:
 //    Calculate Qe*, Cv*, but not nu.
 //    These are needed to calculate other coefficients
 //    and delta temperatures.
 //------------------------------------------------------------------------//

 template<class DS>
 void P13T<DS>::P1Coeffs::calcStarredFields()
 {
     // Get the electron and ion heat capacities.
    
     ccsf CvElec(spMesh);
     matprops.getElectronSpecificHeat(CvElec);

     ccsf CvIon(spMesh);
     matprops.getIonSpecificHeat(CvIon);

    // We need gamma, the electron-ion coupling coefficient.

     ccsf gamma(spMesh);
     matprops.getElectronIonCoupling(gamma);

    // tmpCoeff is a common term to two calculations.
    // Let's just do it once.
    
     ccsf tmpCoeff(spMesh);
     // tmpCoeff = (gamma*dt) / (CvIon + gamma*dt);
     tmpCoeff = gamma / (CvIon/dt + gamma);

     // Shorthand reference names.

     Assert(spCvStar);
     Assert(spQElecStar);
     
     ccsf &CvStar = *spCvStar;
     ccsf &QElecStar = *spQElecStar;

    // CvStar is one of the results, as well as intermediate.
    
     CvStar = CvElec + CvIon * tmpCoeff;
    
     // Calculate QElecStar (Qe*).

     QElecStar = QElectron + (CvIon/dt * (TIon - TElectron) + QIon) * tmpCoeff;

 }

} // end namespace rtt_3T

//---------------------------------------------------------------------------//
//                              end of P1Coeffs.t.cc
//---------------------------------------------------------------------------//
