//----------------------------------*-C++-*----------------------------------//
// timer.hh
// Scott Turner
// 17 March 1998
//---------------------------------------------------------------------------//
// @> Simple timing routine.
//---------------------------------------------------------------------------//

//
// Routine for CPU and Wallclock timing
//

#include "sn/precision.hh"

#include <sys/types.h>
#include <time.h>

void timer( REAL &cpu, REAL &wall )
{

  clock_t cpu_temp;
  time_t  wall_temp;
  time_t  *tloc;

  tloc = 0;

  cpu_temp = clock();
  wall_temp = time( tloc );

  cpu = cpu_temp / 1000000.0;
  wall = wall_temp;

}

//---------------------------------------------------------------------------//
//                              end of timer.hh
//---------------------------------------------------------------------------//

