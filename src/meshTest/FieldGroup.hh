//----------------------------------*-C++-*----------------------------------//
// FieldGroup.hh
// Randy M. Roberts
// Mon Aug 23 17:43:11 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_FieldGroup_hh__
#define __meshTest_FieldGroup_hh__

namespace rtt_meshTest
{
 
//===========================================================================//
// class FieldGroup - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class XD_, class XI_, class XL_, class XDC_>
struct FieldGroup 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef XD_ XD;
    typedef XI_ XI;
    typedef XL_ XL;
    typedef XDC_ XDC;

};

} // end namespace rtt_meshTest

#endif                          // __meshTest_FieldGroup_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/FieldGroup.hh
//---------------------------------------------------------------------------//
