//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRnd_Control_Inline.cc
 * \author Paul Henning
 * \brief  Rnd_Control test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "rng_test.hh"
#include "ds++/Release.hh"
#include "../Random_Inline.hh"
#include "../LFG.h"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>


using rtt_rng::Rnd_Control;
using rtt_rng::LF_Gen;
using rtt_rng::LF_Gen_Ref;

using namespace std;

int seed = 2452423;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void check_basics()
{
    // Test RNG.cc's errprint function
    cerr << "Testing RNG.cc's errprint function "
         << "(a fake message should appear next)." << endl;
    errprint("noway","check_basics","error message goes here");
        
    return;
}

void control_test()
{
    cout << "\nStarting control tests ..." << endl;
    
    // make a controller
    Rnd_Control control(seed);

    // checks
    if (control.get_number() != 1000000000) ITFAILS;
    if (control.get_seed()   != 2452423)    ITFAILS;
    if (control.get_num()    != 0)          ITFAILS;

    // make some random numbers
    LF_Gen r0; control.initialize(r0);
    if (control.get_num()    != 1)          ITFAILS;
    LF_Gen r1; control.initialize(r1);
    if (control.get_num()    != 2)          ITFAILS;
    LF_Gen r2; control.initialize(r2);
    if (control.get_num()    != 3)          ITFAILS;

    LF_Gen rr2; control.initialize(2, rr2);
    if (control.get_num()    != 3)          ITFAILS;

    LF_Gen rr1; control.initialize(1, rr1);
    if (control.get_num()    != 2)          ITFAILS;

    control.set_num(0);

    LF_Gen rr0; control.initialize(rr0);
    if (control.get_num()    != 1)          ITFAILS;

    for (int i = 0; i < 100; i++)
    {
        double rn0  = r0.ran();
        double rrn0 = rr0.ran();
        double rn1  = r1.ran();
        double rrn1 = rr1.ran();
        double rn2  = r2.ran();
        double rrn2 = rr2.ran();
    
        if (rn0 != rrn0)         ITFAILS;
        if (rn1 != rrn1)         ITFAILS;
        if (rn2 != rrn2)         ITFAILS;
    
        if (rn0 == rrn1)         ITFAILS;
        if (rn1 == rrn2)         ITFAILS;
        if (rn2 == rrn0)         ITFAILS;
    }

    if (rtt_rng_test::passed)
    PASSMSG("Rnd_Control simple test ok.");
}

void check_accessors(void)
{
    cout << "\nStarting additional tests...\n" << endl;
    
    // make a controller
    Rnd_Control control(seed);

    // make some random numbers
    LF_Gen r0;
    control.initialize(r0);

    {
        LF_Gen_Ref gr0 = r0.ref();
        double rn0     = gr0.ran();
        if( rn0 < 0 || rn0 >1 )         ITFAILS;
        if( gr0.get_num() != 0 )        ITFAILS;
        if ( LF_Gen::size_bytes() != LFG_DATA_SIZE*sizeof(unsigned int)) ITFAILS;
        LF_Gen sgr0;
        gr0.spawn(sgr0);
        if ( gr0.is_alias_for(sgr0) )                        ITFAILS;
        if ( gr0.get_unique_num() == sgr0.get_unique_num() ) ITFAILS;
    }

    { // test ctors
        
        // create some data
        size_t N(LFG_DATA_SIZE);
        vector<unsigned int> foo(N,0);
        for( size_t i=0; i<N; ++i)
            foo[i] = i+100;

        // some entries in LFG data have extra-special significance; for
        // example, data[OFFSET_hidx] is "hidx".  This must be less than
        // VALID_LM1, which is either 16 or 30, depending on LFG_PARAM_SET.
#if LFG_PARAM_SET == 1
        foo[LFG_DATA_SIZE-4] = 16;
#elif LFG_PARAM_SET == 2
        foo[LFG_DATA_SIZE-4] = 30;
#endif

        // try ctor form 2
        LF_Gen r2( seed, 0 );        

        // try ctor form 3 (foo must have length = LFG_DATA_SIZE
        LF_Gen r3( &foo[0] );

        // try to spawn
        r3.spawn( r2 );
        if( r3 == r2 )
            FAILMSG("LF_Gen spawn creates a distinct LF_Gen object.")
        else
            PASSMSG("LF_Gen equality operator works.")

        // "get_num" returns the generator stream ID.  It is the same 
        // after the call to spawn.
        if (r3.get_num() == r2.get_num()) 
            PASSMSG("LF_Gen spawn creates objects with the same stream IDs.")
        else
            FAILMSG("LF_Gen spawn changed the stream ID number.")

        // Check if the "unique number" is unique enough
        // Note: This number is a combination of the generator 
        // ID and an unsigned int from the state.  With some low 
        // probability it will not be unique, but we only need it 
        // to be unique enough to prevent the resurrection of
        // two particles with the same stream ID, as described in
// https://tf.lanl.gov/sf/go/artf23409?nav=1&_pagenum=1&returnUrlKey=1333470445809
        for (unsigned int i=0; i<1000000; ++i)
        {
            if (r3.get_unique_num() == r2.get_unique_num())
               FAILMSG("LF_Gen unique state number is not unique.")
            r2.ran();
        }

        // Check the id for this stream
        double rn = r3.ran();
        cout << "LF_Gen r3 returns ran() = " << rn << endl;
        unsigned int id = r3.get_num();
        // The value returned by LF_Gen::get_num() lives at the end of the LFG
        // data array.
        if( id != foo[LFG_DATA_SIZE-1] )          ITFAILS;
        if( r3.size() != LFG_DATA_SIZE )          ITFAILS;

        int count(0);
        for( LF_Gen::iterator it=r3.begin(); it != r3.end(); ++it )
            count++;
        if( count != LFG_DATA_SIZE )              ITFAILS;
        
        
        r3.finish_init();
    }
    
    if (rtt_rng_test::passed)
    PASSMSG("Rnd_Control simple test ok.");
    return;
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
        check_basics();
        control_test();
        check_accessors();
    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "While testing tstRnd_Control, " << ass.what()
             << endl;
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_rng_test::passed) 
    {
        cout << "**** tstRnd_Control Test: PASSED" 
             << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstRnd_Control." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstRnd_Control.cc
//---------------------------------------------------------------------------//
