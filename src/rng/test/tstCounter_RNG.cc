//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstCounter_RNG.cc
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Counter_RNG tests.
 * \note   Copyright (C) 2012-2014 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/Counter_RNG.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

#define ITFAILS ut.failure(__LINE__)

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstCounter_RNG(UnitTest &ut)
{
    rtt_rng::Counter_RNG cbrng(12345, 0, 0);

    for (int i = 0; i < 10; ++i)
    {
        std::cout << cbrng.ran() << std::endl;
    }

    ut.passes("Success");
}

void random_existence(UnitTest &ut){
  //test Counter_RNG class
    //create some example data (force it to be uint64_t to test 64-bit)
    vector<uint64_t> foo(4,0);
    foo[0] = 0xcafecafecafecafe;
    foo[1] = 0xbaddecafbaddecaf;
    foo[2] = 0xfabaceaefabaceae;
    foo[3] = 0xdeadbeefdeadbeef;
    //create a Counter_RNG from some data
    Counter_RNG aynRand0(&foo[0]);
    double r0n0 = aynRand0.ran();
    cout << r0n0 << endl;
    //check function get_num()
    if(aynRand0.get_num() != 0xfabaceaefabaceae) ITFAILS;
    //check influence of first argument in data
    foo[0] = 0;
    Counter_RNG aynRand1(&foo[0]);
    double r1n0 = aynRand1.ran();
    cout << r1n0 << endl;
    if(soft_equiv(r1n0, r0n0)) ITFAILS;
    //check influence of second argument in data
    foo[1] = 0;
    Counter_RNG aynRand2(&foo[0]);
    double r2n0 = aynRand2.ran();
    cout << r2n0 << endl;
    if(soft_equiv(r2n0, r1n0)) ITFAILS;
    //check influence of third argument in data
    foo[2] = 0;
    Counter_RNG aynRand3(&foo[0]);
    double r3n0 = aynRand3.ran();
    cout << r3n0 << endl;
    if(soft_equiv(r3n0, r2n0)) ITFAILS;
    //check influence of fourth argument in data
    foo[3] = 0;
    Counter_RNG aynRand4(&foo[0]);
    double r4n0 = aynRand4.ran();
    cout << r4n0 << endl;
    if(soft_equiv(r4n0, r3n0)) ITFAILS;
    //check basic randomness
    double r0n1 = aynRand0.ran();
    double r0n2 = aynRand0.ran();
    cout << r0n1 << endl;
    cout << r0n2 << endl;
    if(soft_equiv(r0n0, r0n1)) ITFAILS;
    if(soft_equiv(r0n0, r0n2)) ITFAILS;
    if(soft_equiv(r0n1, r0n2)) ITFAILS;
    //check alternate (seed,streamnum) constructor
    Counter_RNG aynRand5(0xcafe,0xdecaf, 0);
    if(aynRand5.get_num() != 0xdecaf) ITFAILS;
    //check counter default value
    if(aynRand5.begin()[0] != 0) ITFAILS;
    if(aynRand5.begin()[1] != 0xcafe) ITFAILS;
    //check counter after incrementing
    aynRand5.ran();
    if(aynRand5.begin()[0] != 1) ITFAILS;
    //check counter carry-over
    aynRand5.begin()[0] = 0xffffffffffffffff;
    aynRand5.ran();
    if(aynRand5.begin()[0] != 0) ITFAILS;
    if(aynRand5.begin()[1] != 0xcaff) ITFAILS;
    //check function size()
    if (aynRand0.size() != 4) ITFAILS;
    //check function size_bytes()
    if (aynRand0.size_bytes() != 4 * sizeof(uint64_t)) ITFAILS;
  //test Counter_RNG_Ref class
    //reset the foo values
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //create a Counter_RNG_Ref from the data
    Counter_RNG_Ref aynRandRef(&foo[0],(&foo[0] + 4));
    double rf0n0 = aynRandRef.ran();
    cout << rf0n0 << endl;
    //check function get_num()
    if(aynRand0.get_num() != 0xfabaceaefabaceae) ITFAILS;
    //check function get_unique_num()
    if(aynRand0.get_unique_num() != 0xfabaceaefabaceae) ITFAILS;
    //check influence of first argument in data
    foo[0] = 0;
    double rf0n1 = aynRandRef.ran();
    cout << rf0n1 << endl;
    if(soft_equiv(rf0n1, rf0n0)) ITFAILS;
    //check influence of second argument in data
    foo[0] = 0;
    foo[1] = 0;
    double rf0n2 = aynRandRef.ran();
    cout << rf0n2 << endl;
    if(soft_equiv(rf0n2, rf0n1)) ITFAILS;
    //check influence of third argument in data
    foo[0] = 0;
    foo[2] = 0;
    double rf0n3 = aynRandRef.ran();
    cout << rf0n3 << endl;
    if(soft_equiv(rf0n3, rf0n2)) ITFAILS;
    //check influence of fourth argument in data
    foo[0] = 0;
    foo[3] = 0;
    double rf0n4 = aynRandRef.ran();
    cout << rf0n4 << endl;
    if(soft_equiv(rf0n4, rf0n3)) ITFAILS;
    //reset the foo values
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //check basic randomness
    double rf0n5 = aynRand0.ran();
    double rf0n6 = aynRand0.ran();
    cout << rf0n5 << endl;
    cout << rf0n6 << endl;
    if(soft_equiv(rf0n0, rf0n5)) ITFAILS;
    if(soft_equiv(rf0n0, rf0n6)) ITFAILS;
    if(soft_equiv(rf0n5, rf0n6)) ITFAILS;
    //check counter after incrementing
    foo[0] = 0;
    aynRandRef.ran();
    if(foo[0] != 1) ITFAILS;
    //check counter carry-over
    foo[0] = 0xffffffffffffffff;
    foo[1] = 0;
    aynRandRef.ran();
    if(foo[0] != 0) ITFAILS;
    if(foo[1] != 1) ITFAILS;
  //check Counter_RNG and Counter_RNG_Ref relations
    //reset the foo values
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //check spawning capability
    aynRandRef.spawn(aynRand0);
    aynRand0.spawn(aynRand1);
    if (aynRandRef.get_unique_num() != 0xfabaceae) ITFAILS;
    if (aynRand0.get_unique_num() != 0xfabaceaf) ITFAILS;
    if (aynRand1.get_unique_num() != 0xfabaceb0) ITFAILS;
    if (aynRand0.begin()[0] != 0) ITFAILS;
    if (aynRand0.begin()[2] != 0xfabaceaf) ITFAILS;
    if (aynRand1.begin()[0] != 0) ITFAILS;
    if (aynRand1.begin()[2] != 0xfabaceb0) ITFAILS;
    //check aliasing capability
    Counter_RNG_Ref aynRandAlias0(aynRand0.begin(),aynRand0.end());
    Counter_RNG_Ref aynRandAlias1 = aynRand0.ref();
    if (!aynRandAlias0.is_alias_for(aynRand0)) ITFAILS;
    if (!aynRandAlias1.is_alias_for(aynRand0)) ITFAILS;
    //check equality operators
    Counter_RNG aynRand6(42,42, 0);
    aynRand6 = aynRand0;
    if (!(aynRand0 == aynRand6)) ITFAILS;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        tstCounter_RNG(ut);
        random_existence(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstCounter_RNG, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstCounter_RNG, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstCounter_RNG.cc
//---------------------------------------------------------------------------//
