//----------------------------------*-C++-*----------------------------------//
/*! \file   UnitSystem.hh
 *  \author Kelly Thompson
 *  \brief 
 *  \date   Fri Oct 24 15:07:43 2003
 *  \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_units_UnitSystem_hh
#define rtt_units_UnitSystem_hh

#include <cmath>
#include "ds++/Assert.hh"
#include "UnitSystemType.hh"

namespace rtt_units
{

//===========================================================================//
/*! \class UnitSystem
 *
 * \brief Provide units standardization for Draco. Unit conversion factors
 * defined so that (value in user units) * xxxConversion = (value in SI
 * units) 
 *
 * SI units are:
 *    length  (L) - meters
 *    mass    (M) - kilograms
 *    time    (t) - seconds
 *    temp    (T) - Kelvin
 *    current (I) - Amp
 *    angle (rad) - Radian
 *    quantity    - mole
 *
 * Derived units are:
 *    area             :         L^2
 *    volume           :         L^3
 *    velocity (v)     :         L   / t           
 *    acceleration (a) :         L   / t^2
 *    density (rho)    :     M / L^3
 *    angular velocity : rad     / t
 *    force (f)        :     M * L   / t^2
 *    power (P)        :     M * L^2 / t^3 
 *    torque           :     M * L^2 / t^2
 *    pressure (p)     :     M / L   / t^2
 *    energy           :     M * L^2 / t^2 
 * 
 *    charge                (q) : I * t     (Coulombs = Amp * sec)
 *    electrical potential  (V) : P / I     (Volt = Watt/Amp)
 *    electrical resistance (R) : V / I     (Ohm  = Volt/Amp)
 *    electrical capacitance (C): q / V     (Farads = Coulomb / Volt)
 *    electrical inductance (H) : V * t / I (Henry  = Volt * sec / Amp)
 *
 *    P = V * I
 *    V = I * R
 *    H = I / V / t
 *    C = I * t / V
 *
 *    luminous flux - lumen (lm)
 *    illuminance   - lux (lx)
 *
 * \sa FundUnit
 * \sa PhysicalConstants
 *  
 * \example test/tstUnitSystemType.cc
 * \example tst/tstUnitSystem.cc
 *
 * // Different ways to construct a UnitSystem
 *
 * using rtt_units::UnitSystemType;
 * using rtt_units::UnitSystem;
 * typedef rtt_units::UnitSystemType::
 *
 * UnitSystem uSI( UnitSystemType().SI() );
 * UnitSystem uX4( UnitSystemType().X4() );
 * UnitSystem uCGS( UnitSystemType().L( rtt_units::L_cm  )
 *                                  .M( rtt_units::M_g   )
 *                                  .t( rtt_units::t_s   )
 *                                  .T( rtt_units::T_keV )
 *                                  .I( rtt_units::I_amp )
 *                                  .A( rtt_units::A_rad )
 *                                  .Q( rtt_units::Q_mol ) );
 * UnitSystem uLM_only( UnitSystemType().L( rtt_units::L_cm )
 *                                      .M( rtt_units::M_g  );
 */
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class UnitSystem
{
  public:

    // FRIENDS

//     //! Define the divisor operator for Units/Units.
//     friend UnitSystem operator/( UnitSystem const & op1, 
// 				 UnitSystem const & op2 );

    //! Define the equality operator for Units==Units.
    friend bool operator==( UnitSystem const & op1, 
			    UnitSystem const & op2 );

    //! Define the equality operator for Units==Units.
    friend bool operator!=( UnitSystem const & op1, 
			    UnitSystem const & op2 );
    // CREATORS
    
    //! Prefered constructor
    UnitSystem( UnitSystemType const & ust )
	: d_ust( ust )
    { Require( validUnits() ); }

    //! Default constructor provides SI Units
    UnitSystem()
	: d_ust( UnitSystemType().SI() )
    { Require( validUnits() ); }
    
    // ACCESSORS

    //! Return the conversion factor for length that when multiplied against
    //! values in user-units yields the value in SI-units.
    double      L()     const { return d_ust.L().cf(); }
    std::string Lname() const { return d_ust.L().label(); }

    //! Return the conversion factor for mass that when multiplied against
    //! values in user-units yields the value in SI-units.
    double      M()     const { return d_ust.M().cf(); }
    std::string Mname() const { return d_ust.M().label(); }

    //! Return the conversion factor for time that when multiplied against
    //! values in user-units yields the value in SI-units.
    double t() const { return d_ust.t().cf(); }
    std::string tname() const { return d_ust.t().label(); }

    //! Return the conversion factor for temperature that when multiplied
    //! against values in user-units yields the value in SI-units.
    double      T()     const { return d_ust.T().cf(); }
    std::string Tname() const { return d_ust.T().label(); }

    //! Return the conversion factor for electric current that when
    //! multiplied against values in user-units yields the value in
    //! SI-units. 
    double      I()     const { return d_ust.I().cf(); }
    std::string Iname() const { return d_ust.I().label(); }

    //! Return the conversion factor for angle that when multiplied
    //! against values in user-units yields the value in SI-units.
    double      A()     const { return d_ust.A().cf(); }
    std::string Aname() const { return d_ust.A().label(); }

    //! Return the conversion factor for quantity that when multiplied
    //! against values in user-units yields the value in SI-units.
    double      Q()     const { return d_ust.Q().cf(); }
    std::string Qname() const { return d_ust.Q().label(); }

    //! velocity
    double v() const { return this->L() / this->t(); }
    //! acceleration
    double a() const { return this->L()/std::pow(this->t(),2); }
    //! force
    double f() const { return this->M() * this->L()/std::pow(this->t(),2); }
    //! energy
    double e() const { return this->M() * std::pow(this->L(),2)/pow(this->t(),2); }
    //! power
    double p() const { return this->M() * std::pow(this->L(),2)/pow(this->t(),3); }

    // CLASS IMPLEMENTATION
    
    //! Check whether conversion units are in an acceptable range.
    bool validUnits() const;
    
  private:
    
    // DATA

    // The actual UnitSystem is immutable once created.
    UnitSystemType const d_ust;

}; // end of class UnitSystem

} // end namespace rtt_units

#endif  // rtt__units_UnitSystem_hh

//---------------------------------------------------------------------------//
//                         end of UnitSystem.hh
//---------------------------------------------------------------------------//
