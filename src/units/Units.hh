//----------------------------------*-C++-*----------------------------------//
// Units.hh
// Randy M. Roberts
// Tue Mar 17 14:52:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __units_Units_hh__
#define __units_Units_hh__

#include "ds++/Assert.hh"
#include "units/PhysicalConstants.hh"

#include <cmath>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class Units - 
//
// Date created : Tue Mar 17 14:52:37 1998
// Purpose      : Provide units standardization for Draco.
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

class Units
{
    // FRIENDS

    friend Units operator/(const Units &op1, const Units &op2);

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
    // Convertxxx:
    //    Convert function argument from user units and return it in SI units.
    //------------------------------------------------------------------------//

    template <class FT>
    inline FT ConvertLength(const FT &length) const;
    template <class FT>
    inline FT ConvertMass(const FT &mass) const;
    template <class FT>
    inline FT ConvertTime(const FT &time) const;
    template <class FT>
    inline FT ConvertTemperature(const FT &temperature) const;

    template <class FT>
    inline FT ConvertLength(const FT &length, int n) const;
    template <class FT>
    inline FT ConvertMass(const FT &mass, int n) const;
    template <class FT>
    inline FT ConvertTime(const FT &time, int n) const;
    template <class FT>
    inline FT ConvertTemperature(const FT &temperature, int n) const;
    template <class FT>
    inline FT ConvertTemperature(const FT &temperature, double dn) const;

    template <class FT>
    inline FT ConvertVelocity(const FT &velocity) const;
    template <class FT>
    inline FT ConvertDensity(const FT &density) const;
    template <class FT>
    inline FT ConvertEnergy(const FT &energy) const;

    //------------------------------------------------------------------------//
    // Invertxxx:
    //    Convert function argument from SI and return it in user units.
    //------------------------------------------------------------------------//

    template <class FT>
    inline FT InvertLength(const FT &length) const;
    template <class FT>
    inline FT InvertMass(const FT &mass) const;
    template <class FT>
    inline FT InvertTime(const FT &time) const;
    template <class FT>
    inline FT InvertTemperature(const FT &temperature) const;

    template <class FT>
    inline FT InvertLength(const FT &length, int n) const;
    template <class FT>
    inline FT InvertMass(const FT &mass, int n) const;
    template <class FT>
    inline FT InvertTime(const FT &time, int n) const;
    template <class FT>
    inline FT InvertTemperature(const FT &temperature, int n) const;
    template <class FT>
    inline FT InvertTemperature(const FT &temperature, double dn) const;

    template <class FT>
    inline FT InvertVelocity(const FT &velocity) const;
    template <class FT>
    inline FT InvertDensity(const FT &density) const;
    template <class FT>
    inline FT InvertEnergy(const FT &energy) const;

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
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertLength(const FT &length) const
{
    return length * lengthConversion;
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertMass(const FT &mass) const
{
    return mass * massConversion;
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertTime(const FT &time) const
{
    return time * timeConversion;
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertTemperature(const FT &temperature) const
{
    return temperature * temperatureConversion;
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertLength(const FT &length, int n) const
{
    return length * std::pow(lengthConversion, n);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertMass(const FT &mass, int n) const
{
    return mass * std::pow(massConversion, n);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertTime(const FT &time, int n) const
{
    return time * std::pow(timeConversion, n);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertTemperature(const FT &temperature, int n) const
{
    return temperature * std::pow(temperatureConversion, n);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertTemperature(const FT &temperature, double dn) const
{
    return temperature * std::pow(temperatureConversion, dn);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertLength(const FT &length) const
{
    return length / lengthConversion;
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertMass(const FT &mass) const
{
    return mass / massConversion;
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertTime(const FT &time) const
{
    return time / timeConversion;
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertTemperature(const FT &temperature) const
{
    return temperature / temperatureConversion;
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertLength(const FT &length, int n) const
{
    return length / std::pow(lengthConversion, n);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertMass(const FT &mass, int n) const
{
    return mass / std::pow(massConversion, n);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertTime(const FT &time, int n) const
{
    return time / std::pow(timeConversion, n);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertTemperature(const FT &temperature, int n) const
{
    return temperature / std::pow(temperatureConversion, n);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertTemperature(const FT &temperature, double dn) const
{
    return temperature / std::pow(temperatureConversion, dn);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertVelocity(const FT &velocity) const
{
    return velocity * (lengthConversion / timeConversion);
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertVelocity(const FT &velocity) const
{
    return velocity * (timeConversion / lengthConversion);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertDensity(const FT &density) const
{
    return density * (massConversion /
		      (lengthConversion*lengthConversion*lengthConversion));
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertDensity(const FT &density) const
{
    return density * ((lengthConversion*lengthConversion*lengthConversion) /
		      massConversion);
}

//---------------------------------------------------------------------------//
// Convertxxx:
//    Convert function argument from user units and return it in SI units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::ConvertEnergy(const FT &energy) const
{
    return energy * (massConversion * lengthConversion*lengthConversion /
		     (timeConversion*timeConversion));
}

//---------------------------------------------------------------------------//
// Invertxxx:
//    Convert function argument from SI and return it in user units.
//---------------------------------------------------------------------------//

template <class FT>
inline FT Units::InvertEnergy(const FT &energy) const
{
    return energy * ((timeConversion*timeConversion) /
		     (massConversion * lengthConversion*lengthConversion));	
}

END_NS_XTM  // namespace XTM

#endif                          // __units_Units_hh__

//---------------------------------------------------------------------------//
//                              end of units/Units.hh
//---------------------------------------------------------------------------//
