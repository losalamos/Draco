//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/FieldGroup.hh
 * \author Randy M. Roberts
 * \date   Mon Aug 23 17:43:11 1999
 * \brief  Header file for the FieldGroup class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_FieldGroup_hh__
#define __meshTest_FieldGroup_hh__

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class FieldGroup
 * \brief This class exists because our compiler does not yet support
 *        template template arguments.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class XD_, class XI_, class XL_, class XDC_, int ID>
struct FieldGroup 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef XD_ XD;
    typedef XI_ XI;
    typedef XL_ XL;
    typedef XDC_ XDC;

    static int Id() { return ID; }

};

} // end namespace rtt_meshTest

#endif                          // __meshTest_FieldGroup_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/FieldGroup.hh
//---------------------------------------------------------------------------//
