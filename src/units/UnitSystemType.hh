//----------------------------------*-C++-*----------------------------------//
/*! \file   UnitSystemType.hh
 *  \author Kelly Thompson
 *  \brief  Aggregates a collection of FundUnits to create a complete 
 *          UnitSystemType.
 *  \date   Fri Oct 24 15:04:41 2003
 *  \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __units_UnitSystemType_hh__
#define __units_UnitSystemType_hh__

#include "FundUnit.hh"
 
namespace rtt_units
{

//============================================================================//
/*! 
 * \class UnitSystemType
 *
 * \brief Aggregates a collection of FundUnits to create a complete 
 *        UnitSystemType.
 *
 * \sa UnitSystem
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
 * UnitSystem uCGS( UnitSystemType().L( rtt_units::L_cm )
 *                                  .M( rtt_units::M_g )
 *                                  .t( rtt_units::t_s ) );
 *
 */
//============================================================================//

class UnitSystemType
{
  public:

    // CONSTRUCTORS AND DESTRUCTOR
    
    //! default constructor
    UnitSystemType()
	: d_L( FundUnit< Ltype >( L_null, L_cf, L_labels ) ),
	  d_M( FundUnit< Mtype >( M_null, M_cf, M_labels ) ),
	  d_t( FundUnit< ttype >( t_null, t_cf, t_labels ) ),
	  d_T( FundUnit< Ttype >( T_null, T_cf, T_labels ) ),
	  d_I( FundUnit< Itype >( I_null, I_cf, I_labels ) ),
	  d_A( FundUnit< Atype >( A_null, A_cf, A_labels ) ),
	  d_Q( FundUnit< Qtype >( Q_null, Q_cf, Q_labels ) )
    { /* empty */ }
    
    //! qualified constructor
    UnitSystemType( Ltype myL, Mtype myM, ttype myt, Ttype myT, Itype myI, 
		    Atype myA, Qtype myQ )
	: d_L( FundUnit< Ltype >( myL, L_cf, L_labels ) ),
	  d_M( FundUnit< Mtype >( myM, M_cf, M_labels ) ),
	  d_t( FundUnit< ttype >( myt, t_cf, t_labels ) ),
	  d_T( FundUnit< Ttype >( myT, T_cf, T_labels ) ),
	  d_I( FundUnit< Itype >( myI, I_cf, I_labels ) ),
	  d_A( FundUnit< Atype >( myA, A_cf, A_labels ) ),
	  d_Q( FundUnit< Qtype >( myQ, Q_cf, Q_labels ) )
    { /* empty */ }

    //! Copy constructor
    UnitSystemType( UnitSystemType const & rhs )
	: d_L( rhs.L() ), d_M( rhs.M() ),
	  d_t( rhs.t() ), d_T( rhs.T() ),
	  d_I( rhs.I() ), d_A( rhs.A() ),
	  d_Q( rhs.Q() )	    
    { /* empty */ }
    
    // MANIPULATORS
    
    //! Set SI defaults
    UnitSystemType SI() 
    { return UnitSystemType( L_m,   M_kg, t_s, T_K, I_amp, A_rad, Q_mol ); }

    //! Set X4 defaults
    UnitSystemType X4() 
    { return UnitSystemType( L_cm,  M_g, t_sh, T_keV, I_amp, A_rad, Q_mol ); }

    //! Set cgs defaults
    UnitSystemType CGS()
    { return UnitSystemType( L_cm, M_g, t_s, T_K, I_amp, A_rad, Q_mol ); }

    //! Set a FundUnit type for this UnitSystem

    UnitSystemType& L( Ltype               myType, 
		       double      const * cf     = L_cf, 
		       std::string const & labels = L_labels )
    { 
	this->d_L = FundUnit< Ltype >( myType, cf, labels ); 
	return *this; 
    }
    UnitSystemType& M( Mtype myType, 
		       double const * cf     = M_cf,
		       std::string const & labels = M_labels ) 
    {
	this->d_M = FundUnit< Mtype >( myType, cf, labels );
	return *this; 
    }
    UnitSystemType& t( ttype myType, 
		       double const * cf     = t_cf,
		       std::string const & labels = t_labels ) 
    {
	this->d_t = FundUnit< ttype >( myType, cf, labels ); 
	return *this; 
    }
    UnitSystemType& T( Ttype myType, 
		       double const * cf     = T_cf,
		       std::string const & labels = T_labels ) 
    {
	this->d_T = FundUnit< Ttype >( myType, cf, labels );
	return *this; 
    }
    UnitSystemType& I( Itype myType, 
		       double const * cf     = I_cf,
		       std::string const & labels = I_labels ) 
    { 
	this->d_I = FundUnit< Itype >( myType, cf, labels ); 
	return *this; 
    }
    UnitSystemType& A( Atype myType, 
		       double const * cf     = A_cf,
		       std::string const & labels = A_labels ) 
    {
	this->d_A = FundUnit< Atype >( myType, cf, labels ); 
	return *this; 
    }
    UnitSystemType& Q( Qtype myType, 
		       double const * cf     = Q_cf,
		       std::string const & labels = Q_labels )
    {
	this->d_Q = FundUnit< Qtype >( myType, cf, labels ); 
	return *this; 
    }

    // ACCESSORS

    //! Return a FundUnit type when requested.

    FundUnit< Ltype > L() const { return d_L; }
    FundUnit< Mtype > M() const { return d_M; }
    FundUnit< ttype > t() const { return d_t; }
    FundUnit< Ttype > T() const { return d_T; }
    FundUnit< Itype > I() const { return d_I; }
    FundUnit< Atype > A() const { return d_A; }
    FundUnit< Qtype > Q() const { return d_Q; }
    
  private:

    //! Fundamental unit types.

     FundUnit< Ltype > d_L;
     FundUnit< Mtype > d_M;
     FundUnit< ttype > d_t;
     FundUnit< Ttype > d_T;
     FundUnit< Itype > d_I;
     FundUnit< Atype > d_A;
     FundUnit< Qtype > d_Q;

};

} // end namespace rtt_units

#endif  // __units_UnitSystemType_hh__

//---------------------------------------------------------------------------//
//                     end of UnitSystemType.hh
//---------------------------------------------------------------------------//
