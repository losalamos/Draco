//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstParticle_Stack.cc
 * \author Thomas M. Evans
 * \date   Sat Mar 13 13:49:11 2004
 * \brief  Particle_Stack test.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../Release.hh"
#include "../Particle_Stack.hh"
#include "mc_test.hh"

using namespace std;

using rtt_mc::Particle_Stack;
using rtt_mc::Particle_Containers;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// make a class with an explicit constructor

class DP
{
  private:
    int d_data;
  public:
    explicit DP(int d) : d_data(d) {}
    int data() const { return d_data; }
    void set_data(int d) { d_data = d; }
};

//---------------------------------------------------------------------------//

void stack_test()
{
    // make a particle stack
    Particle_Stack<DP> stack;

    // there should be no particles in it
    if (!stack.empty())    ITFAILS;
    if (stack.size() != 0) ITFAILS;

    // we should catch an illegal pop here
#ifdef REQUIRE_ON
    {
        bool caught = false;
        try
        {
            stack.pop();
        }
        catch (rtt_dsxx::assertion &ass)
        {
            caught = true;
            PASSMSG("Caught illegal pop.");
        }
        if (!caught) ITFAILS;
    }
#endif
    
    // now put some data in
    {
        DP a(1);
        DP b(2);
        DP c(3);
        
        // the stack is last in first out
        stack.push(a);
        stack.push(b);
        stack.push(c);
    }
    
    // check sizes
    if (stack.empty())     ITFAILS;
    if (stack.size() != 3) ITFAILS;

    // check top of stack
    if (stack.top().data() != 3) ITFAILS;
    if (stack.size() != 3)       ITFAILS;

    stack.pop();
    if (stack.size() != 2) ITFAILS;

    // check top of stack
    if (stack.top().data() != 2) ITFAILS;
    if (stack.size() != 2)       ITFAILS;

    stack.pop();
    if (stack.size() != 1) ITFAILS;

    // check top of stack
    if (stack.top().data() != 1) ITFAILS;
    if (stack.size() != 1)       ITFAILS;

    // add another
    stack.push(DP(10));

    if (stack.size() != 2)        ITFAILS;
    if (stack.top().data() != 10) ITFAILS;

    // change the top value
    stack.top().set_data(12);

    // check
    if (stack.size() != 2)        ITFAILS;
    if (stack.top().data() != 12) ITFAILS;
    
    // empty it
    stack.pop();

    // check
    if (stack.size() != 1)       ITFAILS;
    if (stack.top().data() != 1) ITFAILS;

    stack.pop();

    if (!stack.empty())    ITFAILS;
    if (stack.size() != 0) ITFAILS;
    
    if (rtt_mc_test::passed)
        PASSMSG("Particle_Stack ok.");
}

//---------------------------------------------------------------------------//

void subscript_test()
{
    DP a(1);
    DP b(2);
    DP c(3);

    Particle_Stack<DP> stack;
    stack.push(b);
    stack.push(a);
    stack.push(c);

    if (stack.top().data() != 3) ITFAILS;
    if (stack.size() != 3)       ITFAILS;

    // subscript access is first to last
    if (stack[0].data() != 2) ITFAILS;
    if (stack[1].data() != 1) ITFAILS;
    if (stack[2].data() != 3) ITFAILS;

    if (rtt_mc_test::passed)
        PASSMSG("Particle_Stack subscript access ok.");
}

//---------------------------------------------------------------------------//

void containers_test()
{
    SP<DP> a(new DP(1));
    SP<DP> b(new DP(2));
    SP<DP> c(new DP(3));

    // bank and census have smart pointers to the particle type
    Particle_Containers<DP>::Bank   bank;
    Particle_Containers<DP>::Census census;

    bank.push(a);
    bank.push(c);

    census.push(b);

    if (census.size() != 1) ITFAILS;
    if (bank.size() != 2)   ITFAILS;

    if (bank.top()->data() != 3)   ITFAILS;
    if (census.top()->data() != 2) ITFAILS;

    bank.pop();
    bank.pop();
    census.pop();

    if (!census.empty()) ITFAILS;
    if (!bank.empty())   ITFAILS;

    if (rtt_mc_test::passed)
        PASSMSG("Particle_Containers ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " 
		     << rtt_mc::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

        // only run tests on host processor
        if (rtt_c4::node() == 0)
        {
            stack_test();
            subscript_test();
            containers_test();
        }
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tstParticle_Stack, " 
		  << err.what()
		  << std::endl;
	rtt_c4::finalize();
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstParticle_Stack, " 
		  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	if (rtt_mc_test::passed) 
	{
	    std::cout << "**** tstParticle_Stack Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing tstParticle_Stack on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstParticle_Stack.cc
//---------------------------------------------------------------------------//
