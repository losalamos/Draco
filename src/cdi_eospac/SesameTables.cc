//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/SesameTables.cc
 * \author Kelly Thompson
 * \date   Fri Apr  6 08:57:48 2001
 * \brief  Implementation file for SesameTables (mapping material IDs
 *         to Sesame table indexes).
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "SesameTables.hh"
#include "ds++/Assert.hh"

// Need for DEBUG only
#include <iostream>

namespace rtt_cdi_eospac
{

// Constructor.

SesameTables::SesameTables()
    : numReturnTypes( EOS_M_DT ) //  EOS_M_DT = 305 (see eos_Interface.h)
{  
    // Initialize the material map;
    matMap.resize( numReturnTypes, EOS_NullTable );
	    
    // Init a list of return types
    rtMap.resize( numReturnTypes, EOS_NullTable );
}

// Set functions

SesameTables& SesameTables::Pt_DT( unsigned matID )
{
    matMap[ EOS_Pt_DT ] = matID;
    rtMap[ EOS_Pt_DT ] = EOS_Pt_DT;
    return *this;
}

SesameTables& SesameTables::Ut_DT( unsigned matID )
{
    matMap[ EOS_Ut_DT ] = matID;
    rtMap[ EOS_Ut_DT ] = EOS_Ut_DT;
    return *this;
}
SesameTables& SesameTables::T_DPt( unsigned matID ) 
{
    matMap[ EOS_T_DPt ] = matID;
    rtMap[ EOS_T_DPt ] = EOS_T_DPt;
    return *this;
}
SesameTables& SesameTables::T_DUt( unsigned matID ) 
{
    matMap[ EOS_T_DUt ] = matID;
    rtMap[ EOS_T_DUt ] = EOS_T_DUt;
    return *this;
}
SesameTables& SesameTables::Pt_DUt( unsigned matID ) 
{
    matMap[ EOS_Pt_DUt ] = matID;
    rtMap[ EOS_Pt_DUt ] = EOS_Pt_DUt;
    return *this;
}
SesameTables& SesameTables::Ut_DPt( unsigned matID ) 
{
    matMap[ EOS_Ut_DPt ] = matID;
    rtMap[ EOS_Ut_DPt ] = EOS_Ut_DPt;
    return *this;
}
SesameTables& SesameTables::Pic_DT( unsigned matID ) 
{
    matMap[ EOS_Pic_DT ] = matID;
    rtMap[ EOS_Pic_DT ] = EOS_Pic_DT;
    return *this;
}
SesameTables& SesameTables::Uic_DT( unsigned matID ) 
{
    matMap[ EOS_Uic_DT ] = matID;
    rtMap[ EOS_Uic_DT ] = EOS_Uic_DT;
    return *this;
}
SesameTables& SesameTables::T_DPic( unsigned matID ) 
{
    matMap[ EOS_T_DPic ] = matID;
    rtMap[ EOS_T_DPic ] = EOS_T_DPic;
    return *this;
}
SesameTables& SesameTables::T_DUic( unsigned matID ) 
{
    matMap[ EOS_T_DUic ] = matID;
    rtMap[ EOS_T_DUic ] = EOS_T_DUic;
    return *this;
}
SesameTables& SesameTables::Pic_DUic( unsigned matID ) 
{
    matMap[ EOS_Pic_DUic ] = matID;
    rtMap[ EOS_Pic_DUic ] = EOS_Pic_DUic;
    return *this;
}
SesameTables& SesameTables::Uic_DPic( unsigned matID ) 
{
    matMap[ EOS_Uic_DPic ] = matID;
    rtMap[ EOS_Uic_DPic ] = EOS_Uic_DPic;
    return *this;
}
SesameTables& SesameTables::Pe_DT( unsigned matID ) 
{
    matMap[ EOS_Pe_DT ] = matID;
    rtMap[ EOS_Pe_DT ] = EOS_Pe_DT;
    return *this;
}
SesameTables& SesameTables::Ue_DT( unsigned matID ) 
{
    matMap[ EOS_Ue_DT ] = matID;
    rtMap[ EOS_Ue_DT ] = EOS_Ue_DT;
    return *this;
}
SesameTables& SesameTables::T_DPe( unsigned matID ) 
{
    matMap[ EOS_T_DPe ] = matID;
    rtMap[ EOS_T_DPe ] = EOS_T_DPe;
    return *this;
}
SesameTables& SesameTables::T_DUe( unsigned matID ) 
{
    matMap[ EOS_T_DUe ] = matID;
    rtMap[ EOS_T_DUe ] = EOS_T_DUe;
    return *this;
}
SesameTables& SesameTables::Pe_DUe( unsigned matID ) 
{
    matMap[ EOS_Pe_DUe ] = matID;
    rtMap[ EOS_Pe_DUe ] = EOS_Pe_DUe;
    return *this;
}
SesameTables& SesameTables::Ue_DPe( unsigned matID ) 
{
    matMap[ EOS_Ue_DPe ] = matID;
    rtMap[ EOS_Ue_DPe ] = EOS_Ue_DPe;
    return *this;
}
SesameTables& SesameTables::Pc_D( unsigned matID ) 
{
    matMap[ EOS_Pc_D ] = matID;
    rtMap[ EOS_Pc_D ] = EOS_Pc_D;
    return *this;
}
SesameTables& SesameTables::Uc_D( unsigned matID ) 
{
    matMap[ EOS_Uc_D ] = matID;
    rtMap[ EOS_Uc_D ] = EOS_Uc_D;
    return *this;
}
SesameTables& SesameTables::Kr_DT( unsigned matID ) 
{
    matMap[ EOS_Kr_DT ] = matID;
    rtMap[ EOS_Kr_DT ] = EOS_Kr_DT;
    return *this;
}
SesameTables& SesameTables::Keo_DT( unsigned matID )
{
    matMap[ EOS_Keo_DT ] = matID;
    rtMap[ EOS_Keo_DT ] = EOS_Keo_DT;
    return *this;
}
SesameTables& SesameTables::Zfo_DT( unsigned matID )
{
    matMap[ EOS_Zfo_DT ] = matID;
    rtMap[ EOS_Zfo_DT ] = EOS_Zfo_DT;
    return *this;
}
SesameTables& SesameTables::Kp_DT(  unsigned matID )
{
    matMap[ EOS_Kp_DT ] = matID;
    rtMap[ EOS_Kp_DT ] = EOS_Kp_DT;
    return *this;
}
SesameTables& SesameTables::Zfc_DT( unsigned matID )
{
    matMap[ EOS_Zfc_DT ] = matID;
    rtMap[ EOS_Zfc_DT ] = EOS_Zfc_DT;
    return *this;
}
SesameTables& SesameTables::Kec_DT( unsigned matID )
{
    matMap[ EOS_Kec_DT ] = matID;
    rtMap[ EOS_Kec_DT ] = EOS_Kec_DT;
    return *this;
}
SesameTables& SesameTables::Ktc_DT( unsigned matID )
{
    matMap[ EOS_Ktc_DT ] = matID;
    rtMap[ EOS_Ktc_DT ] = EOS_Ktc_DT;
    return *this;
}
SesameTables& SesameTables::B_DT( unsigned matID )
{
    matMap[ EOS_B_DT ] = matID;
    rtMap[ EOS_B_DT ] = EOS_B_DT;
    return *this;
}
SesameTables& SesameTables::Kc_DT( unsigned matID )
{
    matMap[ EOS_Kc_DT ] = matID;
    rtMap[ EOS_Kc_DT ] = EOS_Kc_DT;
    return *this;
}
SesameTables& SesameTables::Tm_D(  unsigned matID )
{
    matMap[ EOS_Tm_D ] = matID;
    rtMap[ EOS_Tm_D ] = EOS_Tm_D;
    return *this;
}
SesameTables& SesameTables::Pm_D(  unsigned matID )
{
    matMap[ EOS_Pm_D ] = matID;
    rtMap[ EOS_Pm_D ] =EOS_Pm_D ;
    return *this;
}
SesameTables& SesameTables::Um_D(  unsigned matID )
{
    matMap[ EOS_Um_D ] = matID;
    rtMap[ EOS_Um_D ] = EOS_Um_D;
    return *this;
}
SesameTables& SesameTables::Tf_D( unsigned matID )
{
    matMap[ EOS_Tf_D ] = matID;
    rtMap[ EOS_Tf_D ] = EOS_Tf_D;
    return *this;
}
SesameTables& SesameTables::Pf_D( unsigned matID )
{
    matMap[ EOS_Pf_D ] = matID;
    rtMap[ EOS_Pf_D ] = EOS_Pf_D;
    return *this;
}
SesameTables& SesameTables::Uf_D( unsigned matID )
{
    matMap[ EOS_Uf_D ] = matID;
    rtMap[ EOS_Uf_D ] = EOS_Uf_D;
    return *this;
}
SesameTables& SesameTables::Gs_D( unsigned matID )
{
    matMap[ EOS_Gs_D ] = matID;
    rtMap[ EOS_Gs_D ] = EOS_Gs_D;
    return *this;
}

//---------------------------------------------------------------------------//
// Get Functions

// Return the enumerated data type associated with the provided integer index
EOS_INTEGER SesameTables::returnTypes( unsigned index ) const
{
    Require( index < numReturnTypes );
    return rtMap[ index ];
}

unsigned SesameTables::matID( EOS_INTEGER returnType ) const
{
    Require( returnType >= 0 );
    return matMap[ returnType ];
}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
// end of SesameTables.cc
//---------------------------------------------------------------------------//
