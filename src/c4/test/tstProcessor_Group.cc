//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstProcessor_Group.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../C4_Functions.hh"
#include "../ParallelUnitTest.hh"
#include "../Release.hh"
#include "../Processor_Group.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
#ifdef C4_MPI
void tstProcessor_Group( UnitTest & ut )
{
    unsigned const pid = rtt_c4::node();
    
    Processor_Group comm(2);
    ut.passes("Processor_Group constructed without exception.");

    unsigned const group_pids = comm.size();
    unsigned const base = pid%2;
    vector<double> sum(1, base+1.);
    comm.sum(sum);
    if (sum[0]==group_pids*(base+1.))
    {
        ut.passes("Correct processor group sum");
    }
    else
    {
        ut.failure("NOT correct processor group sum");
    }
    
    return;
}
#endif // C4_MPI
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    try
    {
        ParallelUnitTest ut( argc, argv, release );
        
#ifdef C4_MPI

       tstProcessor_Group(ut);
#else
       ut.passes("Test inactive for scalar");
#endif //C4_MPI
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                   end of tstProcessor_Group.cc
//---------------------------------------------------------------------------//
