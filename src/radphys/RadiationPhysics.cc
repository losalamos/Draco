//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.cc
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "radphys/RadiationPhysics.hh"

#include "traits/ContainerTraits.hh"
#include "units/PhysicalConstants.hh"

#include <cmath>
#include <limits>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//------------------------------------------------------------------------//
// RadiationPhysics:
//     Constructor based on Units.
//------------------------------------------------------------------------//

RadiationPhysics::RadiationPhysics(const Units &units_)
    : units(units_)
{
    // empty
}

//------------------------------------------------------------------------//
// getStefanBoltzmann:
//     Calculate the Stefan-Boltzmann constant in the user's units.
//------------------------------------------------------------------------//

double RadiationPhysics::getStefanBoltzmann() const
{
    using XTM::PhysicalConstants::stefanBoltzmannSI;
    using std::pow;

    // The Stefan-Boltzmann constant has the following SI units:
    //     W m^-2 K^-4 = kg s^-3 K^-4  (W - watts = kg m^2 s^-3)
    //
    // Therefore, we must convert temperature^-4, time^-3, and
    // mass^1 from SI units to the user's units, using the user's
    // Units.InvertXXX methods.
    
    const double stefanBoltzmann = units.InvertMass(
	units.ConvertTemperature(units.ConvertTime(stefanBoltzmannSI, 3), 4));

    return stefanBoltzmann;
}

//------------------------------------------------------------------------//
// getPlank:
//   Calculate the planck function from the electron temperature.
//   The units for input and output are specified by the units
//   supplied to the class.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getPlanck(const Field &TElectron,
				 Field &planckian) const
{
    Require(ContainerTraits<Field>::conformal(TElectron, planckian));

    const double a = getRadConstant();

    planckian = a * TElectron*TElectron*TElectron*TElectron;
}

//------------------------------------------------------------------------//
// getPlankTemperatureDerivative:
//   Calculate the derivative of the planck function with respect to
//   the electron temperature.
//   The units for input and output are specified by the units
//   supplied to the class.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getPlanckTemperatureDerivative(const Field &TElectron,
						      Field &dplanckdT) const
{
    Require(ContainerTraits<Field>::conformal(TElectron, dplanckdT));

    const double a = getRadConstant();

    dplanckdT = 4.0 * a * TElectron*TElectron*TElectron;
}

//------------------------------------------------------------------------//
// getElectIonCoupling:
//   Calculate the Electron-Ion Coupling function from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectIonCoupling(const Field &density,
					   const Field &TElectron,
					   const Field &z,
					   double abar,
					   Field &electIonCoupling) const
{
    //.... CALCULATES ELECTRON-ION COUPLING COEFFICIENT
    //     USING DEVELOPMENT OF L.D. LANDAU (1936). THIS TREATMENT
    //     ASSUMES THAT THE MATERIAL IS AN IDEAL, CLASSICAL PLASMA, AND
    //     THAT THE ELECTRON AND ION TEMPERATURES ARE WITHIN A
    //     FACTOR OF 1836*abar OF EACH OTHER.

    //.... REFERENCE: LANDAU AND LIFSHITZ, "PHYSICAL KINETICS", 
    //                BUTTERWORTH-HEINEMANN, OXFORD 1995, SECTION 42

    //                MIYAMOTO, "PLASMA PHYSICS FOR NUCLEAR FUSION",
    //                MIT PRESS, 1989. SECTION 4.3

    Require(ContainerTraits<Field>::conformal(density, electIonCoupling));
    Require(ContainerTraits<Field>::conformal(TElectron, electIonCoupling));
    Require(ContainerTraits<Field>::conformal(z, electIonCoupling));

    // eic_constant IS A FUNCTION OF VARIOUS PHYSICAL CONSTANTS IN SI UNITS:
    // eic_constant = (avog_num**3)*sqrt(mass_elect)*(elect_chg**4)*1.e9/ & 
    //                ( sqrt((2.*pi)**3)*(eps_zero**2)*sqrt(boltz_const) )

    // eic_constantSI: (amu**3-m**5-sqrt(K)/kg-s**3)
    
    const double eic_constantSI = 2.993811e+22;

    // Convert eic_constantSI to user's units using InvertXXX().

    const double eic_constant =
	units.InvertLength(
		units.InvertTemperature(
		    units.ConvertMass(units.ConvertTime(eic_constantSI, 3)),
		    0.5),
		5);

    // I must initialize the field with something.
    
    Field lambda_ei = TElectron;

    // Get the Electron-Ion Coulomb Log

    getElectIonCoulombLog(density, TElectron, z, abar, lambda_ei);
    
    electIonCoupling = eic_constant * lambda_ei *
	((z*z*z) / (abar*abar*abar)) * density / pow(TElectron, 1.5);
}

//------------------------------------------------------------------------//
// getElectIonCoulombLog:
//   Calculate the Electron-Ion Coulomb Log (dimensionless) from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectIonCoulombLog(const Field &density,
					     const Field &TElectron,
					     const Field &z, double abar,
					     Field &lambda_ei) const
{
    //... Returns the Coulomb Log (dimensionless) for 
    //    electron-ion interactions. Reference:
    //    NRL Plasma Formulary, 1994. Assumes
    //    (ion_temp) .le. (1836*abar*electron_temp)

    Require(ContainerTraits<Field>::conformal(density, lambda_ei));
    Require(ContainerTraits<Field>::conformal(TElectron, lambda_ei));
    Require(ContainerTraits<Field>::conformal(z, lambda_ei));

    const double Na = PhysicalConstants::avogadro;
    const double eV2K = PhysicalConstants::eV2K;
    
    // Convert from mass density of material to number density of just
    // electrons (in CGS)

    Field neCGS = (Na * 1.e-3 / abar) * z * units.ConvertDensity(density);
    
    // Convert temperature from user units to (Electron-Volts)

    const Field TElect_eV = units.ConvertTemperature(TElectron) / eV2K;

    typedef ContainerTraits<Field> CTF;

    CTF::iterator lit = CTF::begin(lambda_ei);
    CTF::const_iterator Tit = CTF::begin(TElect_eV);
    CTF::const_iterator nit = CTF::begin(neCGS);
    CTF::const_iterator zit = CTF::begin(z);

    while(lit != CTF::end(lambda_ei))
    {
	// Calculate formula switch-over point

	bool lowTemp = *Tit <= 10.0*(*zit)*(*zit);

	if (lowTemp)
	    *lit = 23.0 - log(*zit * sqrt(*nit / pow(*Tit, 3.0)));
	else
	    *lit = 24.0 - log(*zit * sqrt(*nit) / *Tit);

	*lit = (*lit < 1.0) ? 1.0 : *lit;

	lit++; Tit++; nit++; zit++;
    }
}

//------------------------------------------------------------------------//
// getElectronConductionCoeff:
//   Calculate the Electron Conduction Coefficient from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectronConductionCoeff(const Field &density,
						  const Field &TElectron,
						  const Field &z,
						  double abar,
						  Field &electCondCoeff) const
{
    //.... CALCULATES ELECTRON THERMAL CONDUCTION COEFFICIENT
    //     THIS TREATMENT ASSUMES THAT THE MATERIAL IS AN IDEAL, 
    //     CLASSICAL PLASMA.

    //.... REFERENCE: LANDAU AND LIFSHITZ, "PHYSICAL KINETICS", 
    //                BUTTERWORTH-HEINEMANN, OXFORD 1995, SECTION 42

    //                MIYAMOTO, "PLASMA PHYSICS FOR NUCLEAR FUSION",
    //                MIT PRESS, 1989. SECTION 4.3

    Require(ContainerTraits<Field>::conformal(density, electCondCoeff));
    Require(ContainerTraits<Field>::conformal(TElectron, electCondCoeff));
    Require(ContainerTraits<Field>::conformal(z, electCondCoeff));

    //.... ecc_constantSI IS A FUNCTION OF VARIOUS PHYSICAL CONSTANTS
    //     IN SI UNITS:
    //.... ecc_constantSI = (boltz_const**7/2)*sqrt(3)*6*pi*(eps_zero**2) /
    //                      ( sqrt(mass_elect)*(elect_chg**4) )
    //                    
    // ecc_constantSI: 
    
    const double ecc_constantSI = 3.980e-11;

    // Convert ecc_constantSI to user's units using InvertXXX().

    const double ecc_constant =
	units.InvertLength(
	    units.InvertMass(
		units.ConvertTemperature(
		    units.ConvertTime(ecc_constantSI, 3),
		    3.5)));

    // Get the Electron-Electron Coulomb Log

    // I must initialize the field with something.
    
    Field lambda_ee = TElectron;

    getElectElectCoulombLog(density, TElectron, z, abar, lambda_ee);
    
    //  GET ELECTRON GAMMA0

    // I must initialize the field with something.
    
    Field gamma0 = TElectron;

    getElectronGamma0(z, gamma0);

    electCondCoeff = ecc_constant * gamma0* pow(TElectron, 2.5) /
	(lambda_ee * density);
}

//------------------------------------------------------------------------//
// getElectElectCoulombLog:
//   Calculate the Electron-Electron Coulomb Log (dimensionless) from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectElectCoulombLog(const Field &density,
					       const Field &TElectron,
					       const Field &z,
					       double abar,
					       Field &lambda_ee) const
{
    //... Returns the Coulomb Log (dimensionless) for 
    //    electron-electron interactions. Reference:
    //    NRL Plasma Formulary, 1994.

    Require(ContainerTraits<Field>::conformal(density, lambda_ee));
    Require(ContainerTraits<Field>::conformal(TElectron, lambda_ee));
    Require(ContainerTraits<Field>::conformal(z, lambda_ee));

    const double Na = PhysicalConstants::avogadro;
    const double eV2K = PhysicalConstants::eV2K;
    
    // Convert from mass density of material to number density of just
    // electrons (in CGS)

    Field neCGS = (Na * 1.e-3 / abar) * z * units.ConvertDensity(density);
    
    // Convert temperature from user units to (Electron-Volts)

    const Field TElect_eV = units.ConvertTemperature(TElectron) / eV2K;

    typedef ContainerTraits<Field> CTF;

    CTF::iterator lit = CTF::begin(lambda_ee);    
    CTF::const_iterator Tit = CTF::begin(TElect_eV);
    CTF::const_iterator	nit = CTF::begin(neCGS);
    CTF::const_iterator	zit = CTF::begin(z);

    while(lit != CTF::end(lambda_ee))
    {
	// Calculate formula switch-over point

	bool lowTemp = *Tit <= 10.0;

	if (lowTemp)
	    *lit = 23.0 - log(sqrt(*nit / pow(*Tit, 3.0)));
	else
	    *lit = 24.0 - log(sqrt(*nit) / *Tit);

	*lit = (*lit < 1.0) ? 1.0 : *lit;

	lit++; Tit++; nit++; zit++;
    }
}

//------------------------------------------------------------------------//
// getElectronGamma0:
//    Returns a dimensionless constant for use in
//    calculation of electron thermal conduction
//    coefficeint. Correlation developed by
//    C. Cranfill, LANL.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectronGamma0(const Field &z, Field &gamma0) const
{

    Require(ContainerTraits<Field>::conformal(z, gamma0));

    typedef ContainerTraits<Field> CTF;
    
    CTF::const_iterator zit = CTF::begin(z);
    CTF::iterator git = CTF::begin(gamma0);
    
    while (zit != CTF::end(z))
    {
	const double tiny = std::numeric_limits<double>::min();
	
	if ( *zit < tiny )
	{
	    *git = 0.0;
	}
	else
	{
	    using std::sqrt;
	    const double eps = sqrt(2.0) / *zit;
	    const double ke = (433.0 + 180.0*eps)/280.0;
	    const double le = ( 69.0 +  12.0*eps)/(20.0*sqrt(7.0));
	    const double je = ( 13.0 +   4.0*eps)/10.0;

	    // Ratio Of specific heats for ideal monatomic gas
    
	    const double gamma = 5.0/3.0;

	    *git = gamma * ke / ( (gamma-1.0)*(je*ke-le*le) );
	}

	zit++;
	git++;
    }
}

//------------------------------------------------------------------------//
// getIonConductionCoeff:
//   Calculate the Ion Conduction Coefficient from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getIonConductionCoeff(const Field &density,
					     const Field &TElectron,
					     const Field &z,
					     double abar,
					     Field &ionCondCoeff) const
{
    //.... CALCULATES ION THERMAL CONDUCTION COEFFICIENT
    //     THIS TREATMENT ASSUMES THAT THE MATERIAL IS AN IDEAL, 
    //     CLASSICAL PLASMA.

    //.... REFERENCE: LANDAU AND LIFSHITZ, "PHYSICAL KINETICS", 
    //                BUTTERWORTH-HEINEMANN, OXFORD 1995, SECTION 42

    //                MIYAMOTO, "PLASMA PHYSICS FOR NUCLEAR FUSION",
    //                MIT PRESS, 1989. SECTION 4.3

    Require(ContainerTraits<Field>::conformal(density, ionCondCoeff));
    Require(ContainerTraits<Field>::conformal(TElectron, ionCondCoeff));
    Require(ContainerTraits<Field>::conformal(z, ionCondCoeff));

    //.... icc_constantSI IS A FUNCTION OF VARIOUS PHYSICAL CONSTANTS
    //     IN SI UNITS:
    //.... icc_constantSI = (10**1.5)*sqrt(avog_num)*(boltz_const**3.5)
    //                      *sqrt(3)*6*pi*(eps_zero**2)/ (elect_chg**4) 
    //                    
    // icc_constantSI: 
    
    const double icc_constantSI = 9.321e-13;

    // Convert icc_constantSI to user's units using InvertXXX().

    const double icc_constant =
	units.InvertLength(
	    units.InvertMass(
		units.ConvertTemperature(
		    units.ConvertTime(icc_constantSI, 3),
		    3.5)));

    // Get the Ion-Ion Coulomb Log

    // I must initialize the field with something.
    
    Field lambda_ii = TElectron;

    getIonIonCoulombLog(density, TElectron, z, abar, lambda_ii);
    
    //  GET ION GAMMA0

    // I must initialize the field with something.
    
    Field gamma0 = TElectron;

    getIonGamma0(z, TElectron, TElectron, abar, gamma0);

    ionCondCoeff = icc_constant * gamma0* pow(TElectron, 2.5) /
	(lambda_ii * pow(z, 4) * sqrt(abar) * density);
}

//------------------------------------------------------------------------//
// getIonIonCoulombLog:
//   Calculate the Ion-Ion Coulomb Log (dimensionless) from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getIonIonCoulombLog(const Field &density,
					       const Field &TElectron,
					       const Field &z,
					       double abar,
					       Field &lambda_ii) const
{
    //... Returns the Coulomb Log (dimensionless) for 
    //    ion-ion interactions. Reference:
    //    NRL Plasma Formulary, 1994.

    Require(ContainerTraits<Field>::conformal(density, lambda_ii));
    Require(ContainerTraits<Field>::conformal(TElectron, lambda_ii));
    Require(ContainerTraits<Field>::conformal(z, lambda_ii));

    const double Na = PhysicalConstants::avogadro;
    const double eV2K = PhysicalConstants::eV2K;
    
    // Convert from mass density of material to number density of ions (in CGS)

    Field niCGS = (Na * 1.e-3 / abar) * units.ConvertDensity(density);
    
    // Convert temperature from user units to (Electron-Volts)

    const Field TElect_eV = units.ConvertTemperature(TElectron) / eV2K;

    typedef ContainerTraits<Field> CTF;

    CTF::iterator lit = CTF::begin(lambda_ii);
    CTF::const_iterator	Tit = CTF::begin(TElect_eV);
    CTF::const_iterator	nit = CTF::begin(niCGS);
    CTF::const_iterator	zit = CTF::begin(z);

    while(lit != CTF::end(lambda_ii))
    {
	*lit = 23.0 - log(pow(*zit, 3.0)*sqrt(2.0*(*nit) / pow(*Tit, 3.0)));

	*lit = (*lit < 1.0) ? 1.0 : *lit;

	lit++; Tit++; nit++; zit++;
    }
}

//------------------------------------------------------------------------//
// getIonGamma0:
//    Returns a dimensionless constant for use in
//    calculation of ion thermal conduction
//    coefficeint. Correlation developed by
//    C. Cranfill, LANL.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getIonGamma0(const Field &z, const Field &TElect,
				    const Field &TIon, double abar,
				    Field &gamma0) const
{

    Require(ContainerTraits<Field>::conformal(z, gamma0));
    Require(ContainerTraits<Field>::conformal(TElect, gamma0));
    Require(ContainerTraits<Field>::conformal(TIon, gamma0));

    typedef ContainerTraits<Field> CTF;
    
    CTF::const_iterator zit = CTF::begin(z);
    CTF::const_iterator teit = CTF::begin(TElect);
    CTF::const_iterator tiit = CTF::begin(TIon);
    CTF::iterator git = CTF::begin(gamma0);
    
    while (zit != CTF::end(z))
    {
	const double tiny = std::numeric_limits<double>::min();
	
	if ( *zit < tiny || *teit < tiny)
	{
	    *git = 0.0;
	}
	else
	{
	    using std::sqrt;
	    using std::pow;
	    const double eps = pow((*tiit)/( 3672.0*abar*(*teit)), 1.5)
		/ ((*zit)*(*zit));
	    const double ki = (9.0 + 70.0*eps)/14.0;
	    const double li = 3.0/( 5.0*sqrt(7.0) );
	    const double ji = (2.0 +   15.0*eps)/5.0;

	    // Ratio Of specific heats for ideal monatomic gas
    
	    const double gamma = 5.0/3.0;

	    *git = gamma * ki / ( (gamma-1.0)*(ji*ki-li*li) );
	}

	zit++; teit++; tiit++; git++;
    }
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of RadiationPhysics.cc
//---------------------------------------------------------------------------//
