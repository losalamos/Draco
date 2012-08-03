//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstTF_Gen.cc
 * \author Peter Ahrens 
 * \brief  TF_Gen test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "rng_test.hh"
#include "ds++/Release.hh"
#include "TF_Gen.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <vector>

using rtt_rng::TF_Gen;
using rtt_rng::TF_Gen_Ref;
using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void random_existence(){
  //test TF_Gen class
    //create some example data (force it to be uint64_t to test 64-bit)
    vector<unsigned int> foo(4,0);
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //create a TF_Gen from some data
    TF_Gen aynRand0(foo.data());
    double r0n0 = aynRand0.ran();
    cout << r0n0 << endl;
    //check function get_num()
    if(aynRand0.get_num() != 0xdeadbeef) ITFAILS;
    //check influence of first argument in data
    foo[0] = 0;
    TF_Gen aynRand1(foo.data());
    double r1n0 = aynRand1.ran();
    cout << r1n0 << endl;
    if(soft_equiv(r1n0, r0n0)) ITFAILS;
    //check influence of second argument in data
    foo[1] = 0;
    TF_Gen aynRand2(foo.data());
    double r2n0 = aynRand2.ran();
    cout << r2n0 << endl;
    if(soft_equiv(r2n0, r1n0)) ITFAILS;
    //check influence of third argument in data
    foo[2] = 0;
    TF_Gen aynRand3(foo.data());
    double r3n0 = aynRand3.ran();
    cout << r3n0 << endl;
    if(soft_equiv(r3n0, r2n0)) ITFAILS;
    //check influence of fourth argument in data
    foo[3] = 0;
    TF_Gen aynRand4(foo.data());
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
    TF_Gen aynRand5(0xcafebabe,0xbaddecaf);
    if(aynRand5.get_num() != 0xbaddecaf) ITFAILS;
    //check counter default value
    if(aynRand5.begin()[0] != 0) ITFAILS;
    if(aynRand5.begin()[1] != 0) ITFAILS;
    //check counter after incrementing
    aynRand5.ran();
    if(aynRand5.begin()[0] != 1) ITFAILS;
    //check counter carry-over
    aynRand5.begin()[0] = 0xffffffff;
    aynRand5.ran();
    if(aynRand5.begin()[0] != 0) ITFAILS;
    if(aynRand5.begin()[1] != 1) ITFAILS;
    //check function size()
    if (aynRand0.size() != 4) ITFAILS;
    //check function size_bytes()
    if (aynRand0.size_bytes() != 4 * sizeof(unsigned int)) ITFAILS;
  //test TF_Gen_Ref class
    //reset the foo values
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //create a TF_Gen_Ref from the data
    TF_Gen_Ref aynRandRef(foo.data(),(foo.data() + 4));
    double rf0n0 = aynRandRef.ran();
    cout << rf0n0 << endl;
    //check function get_num()
    if(aynRand0.get_num() != 0xdeadbeef) ITFAILS;
    //check function get_unique_num()
    if(aynRand0.get_unique_num() != 0xbaddecafdeadbeef) ITFAILS;
    //check influence of first argument in data
    foo[0] = 0;
    double rf0n1 = aynRandRef.ran();
    cout << rf0n1 << endl;
    if(soft_equiv(rf0n1, rf0n0)) ITFAILS;
    //check influence of second argument in data
    foo[1] = 0;
    double rf0n2 = aynRandRef.ran();
    cout << rf0n2 << endl;
    if(soft_equiv(rf0n2, rf0n1)) ITFAILS;
    //check influence of third argument in data
    foo[2] = 0;
    double rf0n3 = aynRandRef.ran();
    cout << rf0n3 << endl;
    if(soft_equiv(rf0n3, rf0n2)) ITFAILS;
    //check influence of fourth argument in data
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
    foo[0] = 0xffffffff;
    foo[1] = 0;
    aynRandRef.ran();
    if(foo[0] != 0) ITFAILS;
    if(foo[1] != 1) ITFAILS;
  //check TF_Gen and TF_Gen_Ref relations
    //reset the foo values
    foo[0] = 0xcafebabe;
    foo[1] = 0xbaddecaf;
    foo[2] = 0xfabaceae;
    foo[3] = 0xdeadbeef;
    //check spawning capability
    aynRandRef.spawn(aynRand0);
    aynRand0.spawn(aynRand1);
    if (aynRandRef.get_unique_num() != 0xbaddecb0deadbeef) ITFAILS;
    if (aynRand0.get_unique_num() != 0xbaddecb0deadbeef) ITFAILS;
    if (aynRand1.get_unique_num() != 0xbaddecafdeadbeef) ITFAILS;
    if (aynRand0.begin()[0] != 0) ITFAILS;
    if (aynRand0.begin()[2] != 0xfabaceae) ITFAILS;
    if (aynRand1.begin()[0] != 0) ITFAILS;
    if (aynRand1.begin()[2] != 0xfabaceae) ITFAILS;
    //check aliasing capability
    TF_Gen_Ref aynRandAlias0(aynRand0.begin(),aynRand0.end());
    TF_Gen_Ref aynRandAlias1 = aynRand0.ref();
    if (!aynRandAlias0.is_alias_for(aynRand0)) ITFAILS;
    if (!aynRandAlias1.is_alias_for(aynRand0)) ITFAILS;
    //check equality operators
    TF_Gen aynRand6(42,42);
    aynRand6 = aynRand0;
    if (!(aynRand0 == aynRand6)) ITFAILS;
}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
    if (string(argv[arg]) == "--version")
    {
       cout << argv[0] << ": version " << rtt_dsxx::release() 
         << endl;
        return 0;
    }

    cout << "\nThis is rng: version" << rtt_dsxx::release() << "\n" << endl;
    
    try
    {
        // >>> UNIT TESTS
        random_existence();
    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "While testing tstTF_Gen, " << ass.what()
             << endl;
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_rng_test::passed) 
    {
        cout << "**** tstTF_Gen Test: PASSED" 
             << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstTF_Gen." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstTF_Gen.cc
//---------------------------------------------------------------------------//
