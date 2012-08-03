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
#include <vector>
#include <iostream>


using rtt_rng::TF_Gen;
using rtt_rng::TF_Gen_Ref;

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void random_existence(){
    std::vector<unsigned int> foo(4);
    foo[0] = 0xbadcafe0;
    foo[1] = 0xbadcafe1;
    foo[2] = 0xbadcafe2;
    foo[3] = 0xdeadbeef;
    rtt_rng::TF_Gen_Ref aynRandRef(foo.data(),(foo.data() + 4));
    rtt_rng::TF_Gen aynRand0(foo.data());
    rtt_rng::TF_Gen aynRand1(0,0);
    rtt_rng::TF_Gen aynRand2(0,1);
    rtt_rng::TF_Gen aynRand3(42,42);
    aynRandRef.spawn(aynRand1);
    aynRand1.spawn(aynRand2);
    rtt_rng::TF_Gen_Ref aynRandAlias0(aynRand0.begin(),aynRand0.end());
    rtt_rng::TF_Gen_Ref aynRandAlias1 = aynRand0.ref();
    rtt_rng::TF_Gen aynRand4(0,0);
    aynRand4 = aynRand0;
    if (!(aynRand0 == aynRand4)) ITFAILS;

    double randNum0 = aynRandRef.ran();
    double randNum1 = aynRandRef.ran();
    double randNum2 = aynRand1.ran();
    double randNum3 = aynRand2.ran();
    std::cout << randNum0 << std::endl;
    std::cout << randNum1 << std::endl;
    std::cout << randNum2 << std::endl;
    std::cout << randNum3 << std::endl;
    if (randNum0 == 0.0) ITFAILS;
    if (randNum0 == randNum1) ITFAILS;
    if (aynRandRef.get_num() != 0xdeadbeefbadcafe2) ITFAILS;

    if (randNum0 == randNum2) ITFAILS;
    if (aynRand1.get_num() != 0xdeadbef0badcafe2) ITFAILS;
    if (randNum0 == randNum3) ITFAILS;
    if (randNum2 == randNum3) ITFAILS;
    if (aynRand2.get_num() != 0xdeadbef1badcafe2) ITFAILS;
    
    randNum0 = aynRand0.ran();
    randNum1 = aynRand0.ran();
    std::cout << randNum0 << std::endl;
    std::cout << randNum1 << std::endl;
    if (randNum0 == 0.0) ITFAILS;
    if (randNum0 == randNum1) ITFAILS;
    if (aynRand0.get_num() != 0xdeadbeefbadcafe2) ITFAILS;
    if (aynRand0.size() != 4) ITFAILS;
    if (aynRand0.size_bytes() != 4 * sizeof(unsigned int)) ITFAILS;
    if (!aynRandAlias0.is_alias_for(aynRand0)) ITFAILS;
    if (!aynRandAlias1.is_alias_for(aynRand0)) ITFAILS;

    if (aynRand3.get_num() != 0x0000002a0000002a) ITFAILS;
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
