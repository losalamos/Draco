//----------------------------------*-C++-*----------------------------------//
// Units.hh
// Randy M. Roberts
// Tue Mar 17 14:52:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_Units_hh__
#define __3T_Units_hh__

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class Units - 
//
// Date created : Tue Mar 17 14:52:37 1998
// Purpose      : Provide units standardization for the 3T package.
//                Unit conversion factors defined so that
//                (value in user units) * xxxConversion = (value in SI units)
//
//                SI units are:
//                    length - meters
//                    mass   - kilograms
//                    time   - seconds
//                    temp   - Kelvin
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "ds++/Assert.hh"
#include "3T/PhysicalConstants.hh"

class Units
{

    // DATA

    double lengthConversion;
    double massConversion;
    double timeConversion;
    double temperatureConversion;
    
  public:

    // CLASS METHODS

    //------------------------------------------------------------------------//
    // getSIUnits:
    //     A convenience utility to return an SI-based Units object.
    //------------------------------------------------------------------------//
    
    static Units getSIUnits() { return Units(1.0, 1.0, 1.0, 1.0); }
    
    //------------------------------------------------------------------------//
    // getAstroPhysUnits:
    //     A convenience utility to return an Units object
    //     based on AstroPhysical units.
    //------------------------------------------------------------------------//
    
    static Units getAstroPhysUnits()
    {
	// Convenience utility to return Astro-Physical units
	// length - centimeters
	// mass   - grams
	// time   - 10^-8 seconds (1 shake)
	// temp   - keV

	using XTM::PhysicalConstants::eV2K;

	const double mPcm = 1.0e-2;
	const double kgPg = 1.0e-3;
	const double sPsh = 1.0e-8;
	const double KPkeV = 1.0e3*eV2K;
	
	return Units(mPcm, kgPg, sPsh, KPkeV);
    }
    
    // CREATORS
    
    Units(double lengthConversion_, double massConversion_,
	  double timeConversion_, double temperatureConversion_)
	: lengthConversion(lengthConversion_),
	  massConversion(massConversion_),
	  timeConversion(timeConversion_),
	  temperatureConversion(temperatureConversion_)
    {
	Require(validUnits());
    }

    // The default Units constructor yields
    // user-units = SI-units
    
    Units()
	: lengthConversion(1.0),
	  massConversion(1.0),
	  timeConversion(1.0),
	  temperatureConversion(1.0)
    {
	Require(validUnits());
    }
    

    // MANIPULATORS
    // none
    
    // ACCESSORS

    // These getxxxConversion() accessors return the conversion factor that
    // when multiplied against values in user-units yields
    // the value in SI-units.
    
    double getLengthConversion() const { return lengthConversion; }
    double getMassConversion() const { return massConversion; }
    double getTimeConversion() const { return timeConversion; }
    double getTemperatureConversion() const { return temperatureConversion; }

    //------------------------------------------------------------------------//
    // User2SIxxx:
    //    Convert function argument from user units and return it in SI units.
    //------------------------------------------------------------------------//

    inline double User2SILength(double length) const;
    inline double User2SIMass(double mass) const;
    inline double User2SITime(double time) const;
    inline double User2SITemperature(double temperature) const;

    inline double User2SIVelocity(double velocity) const;
    inline double User2SIDensity(double density) const;

    //------------------------------------------------------------------------//
    // SI2Userxxx:
    //    Convert function argument from SI and return it in user units.
    //------------------------------------------------------------------------//

    inline double SI2UserLength(double length) const;
    inline double SI2UserMass(double mass) const;
    inline double SI2UserTime(double time) const;
    inline double SI2UserTemperature(double temperature) const;

    inline double SI2UserVelocity(double velocity) const;
    inline double SI2UserDensity(double density) const;

  protected:
    
    // IMPLEMENTATION

    //------------------------------------------------------------------------//
    // validUnits:
    //    Check whether conversion units are in an acceptable range.
    //------------------------------------------------------------------------//

    inline bool validUnits() const;

    // CLASS IMPLEMENTATION

    //------------------------------------------------------------------------//
    // minConversion:
    //    A class utility that returns the minimum allowed conversion factor.
    //------------------------------------------------------------------------//

    static double minConversion();
};

// INLINE DEFINITIONS

//---------------------------------------------------------------------------//
// validUnits:
//    Check whether conversion units are in an acceptable range.
//---------------------------------------------------------------------------//

inline bool Units::validUnits() const
{
    return lengthConversion >= minConversion()
	&& massConversion >= minConversion()
	&& timeConversion >= minConversion()
	&& temperatureConversion >= minConversion();
}

//---------------------------------------------------------------------------//
// User2SIxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

inline double Units::User2SILength(double length) const
{
    return length * lengthConversion;
}

//---------------------------------------------------------------------------//
// User2SIxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

inline double Units::User2SIMass(double mass) const
{
    return mass * massConversion;
}

//---------------------------------------------------------------------------//
// User2SIxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

inline double Units::User2SITime(double time) const
{
    return time * timeConversion;
}

//---------------------------------------------------------------------------//
// User2SIxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

inline double Units::User2SITemperature(double temperature) const
{
    return temperature * temperatureConversion;
}

//---------------------------------------------------------------------------//
// SI2Userxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

inline double Units::SI2UserLength(double length) const
{
    return length / lengthConversion;
}

//---------------------------------------------------------------------------//
// SI2Userxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

inline double Units::SI2UserMass(double mass) const
{
    return mass / massConversion;
}

//---------------------------------------------------------------------------//
// SI2Userxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

inline double Units::SI2UserTime(double time) const
{
    return time / timeConversion;
}

//---------------------------------------------------------------------------//
// SI2Userxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

inline double Units::SI2UserTemperature(double temperature) const
{
    return temperature / temperatureConversion;
}

//---------------------------------------------------------------------------//
// User2SIxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

inline double Units::User2SIVelocity(double velocity) const
{
    return velocity * lengthConversion / timeConversion;
}

//---------------------------------------------------------------------------//
// SI2Userxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

inline double Units::SI2UserVelocity(double velocity) const
{
    return velocity * timeConversion / lengthConversion;
}

END_NS_XTM  // namespace XTM

#endif                          // __3T_Units_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Units.hh
//---------------------------------------------------------------------------//
