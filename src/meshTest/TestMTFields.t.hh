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
void TestMTFields<MTFactory>::run()
{
    
    // Get the mesh and field constructor from the mesh factory.

    std::string name = Name(FGRP::Id());

    os_m << "TestMTFields::run() running " << name << " test." << std::endl;
    
    MTFactoryProduct meshProduct = meshFactory_m.create();
    MT &mesh = meshProduct.mesh();
    FieldConstructor &fCtor = meshProduct.fieldConstructor();

    FieldTester<MT, FGRP> tester(name, fCtor, mesh, os_m);

    tester.run();
    
    passed_m = tester.passed();
}
 
} // end namespace rtt_meshTest


//---------------------------------------------------------------------------//
//                              end of TestMTFields.t.hh
//---------------------------------------------------------------------------//
