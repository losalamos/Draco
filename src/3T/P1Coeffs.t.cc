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
                              const DiffusionSolver &solver_,
			      const MaterialProperties &matprops_,
			      const ncvsf &velocity_,
			      const RadiationStateField &prevStateField_,
			      const ccsf &QRad_,
			      const ccsf &QElectron_,
			      const ccsf &QIon_,
			      const ccsf &TElectron_,
			      const ccsf &TIon_)
     : p13T(p13T_), spMesh(p13T_.getMesh()),
       dt(dt_), groupNo(groupNo_), options(options_), solver(solver_),
       matprops(matprops_), velocity(velocity_),
       prevStateField(prevStateField_),
       QRad(QRad_), QElectron(QElectron_), QIon(QIon_),
       TElectron(TElectron_), TIon(TIon_),
       xQEEM(spMesh), xREEM(spMesh), xQRadBar(spMesh),
       xQElecStar(spMesh), xCvStar(spMesh), xNu(spMesh),
       xSigmaAbsBar(spMesh), xD(spMesh), xFprime(spMesh),
       xXiTilde(spMesh), xXiBracket(spMesh), xXiOmegaBracket(spMesh)
 {
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
     fcdsf fcSigmaAbs(spMesh);
     matprops.getSigmaAbsorption(groupNo, fcSigmaAbs);

     // Allocate and calculate the diffusion constant.

     xD = (1.0/3.0) / (sigmaTotal + tauP1);

     // Make shorthand references.

     ccsf &QEEM = xQEEM;
     ccsf &REEM = xREEM;
     ccsf &QRadBar = xQRadBar;
     const ccsf &QElecStar = xQElecStar;
     const ccsf &CvStar = xCvStar;
     const ccsf &nu = xNu;
     const ccsf &xiBracket = xXiBracket;
     const fcdsf &xiOmegaBracket = xXiOmegaBracket;

     if (options.getIsCoupledMaterial())
     {
       	 calcVelocityCorrections();
	 calcStarredFieldsAndNu();

	 cerr << "QElecStar: " << *QElecStar.begin() << endl;
	 cerr << "CvStar: " << *CvStar.begin() << endl;
	 cerr << "nu: " << *nu.begin() << endl;
    
	 // Calculate modified sigma absorption

	 xSigmaAbsBar = (1.0 - nu) * sigmaAbs + tau;

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

	 QRadBar = tau*prevStateField.phi + xiBracket + QEEM + QRad;
     }
     else
     {
	 xSigmaAbsBar = sigmaAbs + tau;
	 REEM = 0.0;
	 QEEM = 0.0;
	 QRadBar = tau*prevStateField.phi + QRad;
     }

     // Calculate the "telegraph" term to the P1 equation.

     xFprime = (tauP1*prevStateField.F + xiOmegaBracket) / (sigmaTotal + tauP1);

 }

 //------------------------------------------------------------------------//
 // calcVelocityCorrections:
 //    Calculate the velocity correction terms.
 //------------------------------------------------------------------------//

 template<class DS>
 void P13T<DS>::P1Coeffs::calcVelocityCorrections()
 {
     // Set the radiation physics to the given units.
    
     const RadiationPhysics radPhys(matprops.getUnits());

     // obtain the speed of light

     double c = radPhys.getLightSpeed();

     // get cross sections

     fcdsf fcSigmaAbs(spMesh);
     matprops.getSigmaAbsorption(groupNo, fcSigmaAbs);
     
     fcdsf fcSigmaScatter(spMesh);
     matprops.getSigmaScattering(groupNo, fcSigmaScatter);
     
     // obtain the vertex centered cross sections here

     DiscKineticEnergyField vcSigmaTotal(spMesh);
     DiscKineticEnergyField vcSigmaAbs(spMesh);
     DiscKineticEnergyField vcSigmaScattering(spMesh);

     //fill vcSigma's here

     // a nested scope
     {
	 fcdsf fcSigmaTotal(spMesh);
	 matprops.getSigmaTotal(groupNo, fcSigmaTotal);

	 p13T.mapCrossSections(vcSigmaTotal, fcSigmaTotal);
     }
     p13T.mapCrossSections(vcSigmaAbs, fcSigmaAbs);
     p13T.mapCrossSections(vcSigmaScattering, fcSigmaScatter);
     
     // Obtain unit vectors

     DiscMomentumField::value_type e1;
     DiscMomentumField::value_type e2;
     DiscMomentumField::value_type e3;
     e1(0) = 1.;
     e1(1) = 0.;
     e1(2) = 0.;
     e2(0) = 0.;
     e2(1) = 1.;
     e2(2) = 0.;
     e3(0) = 0.;
     e3(1) = 0.;
     e3(2) = 1.;

     // put the velocities and fluxes into the right field type

     DiscMomentumField vvelocity(spMesh);
     MT::gather ( vvelocity, velocity, MT::OpAssign() );
     DiscKineticEnergyField phi(spMesh);
     MT::gather ( phi, prevStateField.phi, MT::OpAssign() );

     // obtain face normal components

     fcdvsf face_normals( spMesh );
     face_normals = spMesh->get_fn();
     DiscFluxField n1( spMesh );
     DiscFluxField n2( spMesh );
     DiscFluxField n3( spMesh );
     typedef typename fcdvsf::value_type vec;
     {
         DiscFluxField::iterator iter1 = n1.begin();
         DiscFluxField::iterator iter2 = n2.begin();
         DiscFluxField::iterator iter3 = n3.begin();
         for (fcdvsf::iterator iter = face_normals.begin();
              iter != face_normals.end(); iter++)
         {
             *iter1++ = vec::dot(*iter, e1);
             *iter2++ = vec::dot(*iter, e2);
             *iter3++ = vec::dot(*iter, e3);
         }
     }

     // calculate dot products and components of velocity term

     DiscKineticEnergyField velocity_sqrd(spMesh);
     {
         DiscKineticEnergyField::iterator itersq = velocity_sqrd.begin();
         for (DiscMomentumField::iterator iter = vvelocity.begin();
              iter != vvelocity.end(); iter++)
         {
             *itersq++ = vec::dot(*iter, *iter);
         }
     }
     DiscKineticEnergyField velocity1(spMesh);
     DiscKineticEnergyField velocity2(spMesh);
     DiscKineticEnergyField velocity3(spMesh);
     {
         DiscKineticEnergyField::iterator iter1 = velocity1.begin();
         DiscKineticEnergyField::iterator iter2 = velocity2.begin();
         DiscKineticEnergyField::iterator iter3 = velocity3.begin();
         for (DiscMomentumField::iterator iter = vvelocity.begin();
              iter != vvelocity.end(); iter++)
         {
             *iter1++ = vec::dot(*iter, e1);
             *iter2++ = vec::dot(*iter, e2);
             *iter3++ = vec::dot(*iter, e3);
         }
     }

     // determine the vertex to cell volume ratios

     ccsf cell_volumes(spMesh);
     spMesh->get_cell_volumes(cell_volumes);
     DiscKineticEnergyField vertex_volumes(spMesh);
     spMesh->get_vertex_volumes(vertex_volumes);
     DiscKineticEnergyField vc_cell_volumes(spMesh);
     MT::gather ( vc_cell_volumes, cell_volumes, MT::OpAssign() );
     DiscKineticEnergyField vc_volume_ratios(spMesh);
     vc_volume_ratios = vertex_volumes/vc_cell_volumes;

     // Make shorthand references.

     ccsf &xiTilde = xXiTilde;
     ccsf &xiBracket = xXiBracket;
     DiscFluxField &xiOmegaBracket = xXiOmegaBracket;

     // calculate xiTilde

     DiscFluxField sigmaF(spMesh);
     sigmaF = fcSigmaAbs*prevStateField.F;
     DiscKineticEnergyField DKEField(spMesh);
     solver.dotProduct(DKEField, sigmaF, vvelocity);
     DKEField -= (4./(3.*c))*vcSigmaAbs*phi*velocity_sqrd;
     DKEField *= vc_volume_ratios;
     MT::scatter ( xiTilde, DKEField, MT::OpAddAssign() );
     xiTilde *= (2./c);

     // calculate xiBracket

     sigmaF = (fcSigmaAbs - fcSigmaScatter)*prevStateField.F;
     solver.dotProduct(DKEField, sigmaF, vvelocity);
     DKEField -= (4./(3.*c))*(vcSigmaAbs-vcSigmaScattering)*phi*velocity_sqrd;
     DKEField *= vc_volume_ratios;
     MT::scatter ( xiBracket, DKEField, MT::OpAddAssign() );
     xiBracket /= c;

     // calculate xiOmegaBracket

     DiscKineticEnergyField DKEField1(spMesh);
     DiscKineticEnergyField DKEField2(spMesh);
     DiscKineticEnergyField DKEField3(spMesh);
     // note: for orthogonal mesh Vvf/Vf = 1/4
     DKEField1 = vcSigmaTotal*phi*velocity1/4.;
     DKEField2 = vcSigmaTotal*phi*velocity2/4.;
     DKEField3 = vcSigmaTotal*phi*velocity3/4.;
     DiscFluxField DFField1(spMesh);
     DiscFluxField DFField2(spMesh);
     DiscFluxField DFField3(spMesh);
     MT::scatter( DFField1, DKEField1, MT::OpAddAssign() );
     MT::scatter( DFField2, DKEField2, MT::OpAddAssign() );
     MT::scatter( DFField3, DKEField3, MT::OpAddAssign() );
     xiOmegaBracket = (4./(3.*c))*(DFField1*n1 + DFField2*n2 + DFField3*n3);
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

     const ccsf &CvStar = xCvStar;
     ccsf &nu = xNu;

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

     ccsf &CvStar = xCvStar;
     ccsf &QElecStar = xQElecStar;
     const ccsf &xiTilde = xXiTilde;

    // CvStar is one of the results, as well as intermediate.
    
     CvStar = CvElec + CvIon * tmpCoeff;
    
     // Calculate QElecStar (Qe*).

     QElecStar = QElectron + (CvIon/dt * (TIon - TElectron) + QIon) * tmpCoeff
                 - xiTilde;

 }

} // end namespace rtt_3T

//---------------------------------------------------------------------------//
//                              end of P1Coeffs.t.cc
//---------------------------------------------------------------------------//
