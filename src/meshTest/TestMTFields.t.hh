//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTFields.t.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Mon Aug 23 16:49:16 1999
 * \brief  Implementation file for the TestMTFields class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMTFields.hh"
#include "FieldTester.hh"

#include <iostream>

namespace rtt_meshTest
{

template<class MTFactory>
template<class FGRP>
void TestMTFields<MTFactory>::run(const std::string &name_in)
{
    
    // Get the mesh and field constructor from the mesh factory.
    
    MTFactoryProduct meshProduct = meshFactory_m.create();
    MT &mesh = meshProduct.mesh();
    FieldConstructor &fCtor = meshProduct.fieldConstructor();

    FieldTester<MT, FGRP> tester(name_in, fCtor, mesh, os_m);

    tester.run();
    
    passed_m = tester.passed();
}
 
} // end namespace rtt_meshTest


//---------------------------------------------------------------------------//
//                              end of TestMTFields.t.hh
//---------------------------------------------------------------------------//
