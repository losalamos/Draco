//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Serial.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:06:25 2002
 * \brief  Implementation of C4 serial option.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_SCALAR

#include "C4_Serial.hh"
#include <unistd.h>
#include <sys/times.h>
#include <cstdlib>

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// Null source/destination rank
//---------------------------------------------------------------------------//

const int proc_null = -1;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

void initialize(int &argc, char **&argv)
{
}

//---------------------------------------------------------------------------//

void finalize()
{
}

//---------------------------------------------------------------------------//

void free_inherited_comm()
{
}

//---------------------------------------------------------------------------//
// QUERY FUNCTIONS
//---------------------------------------------------------------------------//

int node()
{
    return 0;
}

//---------------------------------------------------------------------------//

int nodes()
{
    return 1;
}

//---------------------------------------------------------------------------//
// BARRIER FUNCTIONS
//---------------------------------------------------------------------------//

void global_barrier()
{
}

//---------------------------------------------------------------------------//
// TIMING FUNCTIONS
//---------------------------------------------------------------------------//

double wall_clock_time()
{
    tms now;
    return times( &now ) / wall_clock_resolution();
}
double wall_clock_time( tms & now )
{
    return times( &now ) / wall_clock_resolution();
}

//---------------------------------------------------------------------------//

double wall_clock_resolution()
{
    return static_cast<double>(sysconf(_SC_CLK_TCK));
}

//---------------------------------------------------------------------------//
// PROBE/WAIT FUNCTIONS
//---------------------------------------------------------------------------//

bool probe(int  /* source */, 
	   int  /* tag */,
	   int &/* message_size */)
{    
    return false;
}

void blocking_probe(int  /* source */, 
                    int  /* tag */,
                    int &/* message_size */)
{    
    Insist(false, "no messages expected in serial programs!");
}

//---------------------------------------------------------------------------//
// ABORT
//---------------------------------------------------------------------------//

int abort(int error)
{
    // This test is not recorded as tested by BullseyeCoverage because abort
    // terminates the execution and BullseyeCoverage only reports coverage for
    // function that return control to main().
    
    // call system exit
    std::abort();
    return error;
}

//---------------------------------------------------------------------------//

bool isScalar()
{
    return true;
}

} // end namespace rtt_c4

#endif // C4_SCALAR

//---------------------------------------------------------------------------//
//                              end of C4_Serial.cc
//---------------------------------------------------------------------------//
