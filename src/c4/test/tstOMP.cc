//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/phw.cc
 * \author Kelly Thompson
 * \date   Tue Jun  6 15:03:08 2006
 * \brief  Parallel application used by the unit test for tstApplicationUnitTest.
 * \note   Copyright (C) 2006-2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../global.hh"
#include "../ParallelUnitTest.hh"
#include "../Timer.hh"
#include "../gatherv.hh"

#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <complex>
#include <cstdio>
#include <vector>
#include <string>

#include "c4_omp.h"

using namespace rtt_c4;

typedef std::complex<double> complex;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#define PASSMSG(A) ut.passes(A)
#define FAILMSG(A) ut.failure(A)
#define ITFAILS    ut.failure( __LINE__ )

//---------------------------------------------------------------------------//

bool topology_report(void)
{
    size_t const mpi_ranks = rtt_c4::nodes();
    size_t const my_mpi_rank = rtt_c4::node();

    // Store proc name on local proc
    std::string my_pname = rtt_c4::get_processor_name();
    Remember( size_t namelen = my_pname.size(); );

    // Create a container on IO proc to hold names of all nodes.    
    std::vector< std::string > procnames(mpi_ranks);
    
    // Gather names into pnames on IO proc.
    rtt_c4::indeterminate_gatherv( my_pname, procnames );

    // Is there only 1 MPI rank per machine node?
    int one_mpi_rank_per_node(0);
    
    // Look at the data found on the IO proc.
    if( my_mpi_rank == 0 )
    {
        Check( procnames[my_mpi_rank].size() == namelen);

        // Count unique processors
        std::vector< std::string > unique_processor_names;
        for( size_t i=0; i<mpi_ranks; ++i )
        {
            bool found(false);
            for( size_t j=0; j<unique_processor_names.size(); ++j )
                if( procnames[i] == unique_processor_names[j] )
                    found = true;
            if( ! found )
                unique_processor_names.push_back( procnames[i] );
        }

        // Print a report
        std::cout << "\nWe are using " << mpi_ranks << " mpi rank(s) on "
                  << unique_processor_names.size() << " unique nodes.";
        
        for( size_t i=0; i<mpi_ranks; ++i )
            std::cout << "\n  - MPI rank " << i <<" is on " << procnames[i];
        std::cout << "\n" << std::endl;

        if( mpi_ranks == unique_processor_names.size() )
            one_mpi_rank_per_node = 1;
    }

    rtt_c4::broadcast( &one_mpi_rank_per_node, 1, 0 );
    
    // return 't' if 1 MPI rank per machine node.
    return (one_mpi_rank_per_node == 1);
}

//---------------------------------------------------------------------------//

void topo_report(rtt_dsxx::UnitTest &ut, bool & one_mpi_rank_per_node )
{
    // Determine if MPI ranks are on unique machine nodes:
    //
    // If there are multiple MPI ranks per machine node, then don't use OMP
    // because OMP can't restrict its threads to running only on an MPI rank's
    // cores.  The OMP threads will be distributed over the whole machine
    // node.  For example, we might choose to use 4 MPI ranks on a machine
    // node with 16 cores.  Ideally, we could allow each MPI rank to use 4 OMP
    // threads for a maximum of 4x4=16 OMP threads on the 16 core node.
    // However, because OMP doesn't know about the MPI ranks sharing the 16
    // cores, the even distribution of OMP threads is not guaranteed.
    //
    // So - if we have more than one MPI rank per machine node, then turn off
    // OMP threads.
    one_mpi_rank_per_node = topology_report();

    std::string procname = rtt_c4::get_processor_name();

#ifdef USE_OPENMP
    int tid(-1);
    int nthreads(-1);

    if( one_mpi_rank_per_node )
        nthreads = omp_get_max_threads();
    else
    {
        // More than 1 MPI rank per node --> turn off OMP.
        nthreads = 1;
        omp_set_num_threads( nthreads );
    }
        
#pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        if( tid == nthreads-1 )
        {
            std::cout << "Using OMP threads."
                      << "\n   MPI node       : " << node()
                      << "\n   MPI max nodes  : " << nodes()
                      << "\n   OMP thread     : " << tid
                      << "\n   OMP max threads: " << nthreads
                      << "\n   procname(IO)   : " << procname
                      << "\n" << std::endl;
        }
        if( tid < 0 || tid >= nthreads )  ITFAILS;
    }
#else
    { // not OMP
        std::cout << "OMP thread use is disabled."
                  << "\n   MPI node       : " << node()
                  << "\n   MPI max nodes  : " << nodes()
                  << "\n   procname(IO)   : " << procname
                  << "\n" << std::endl;
    }
#endif
    
    return;
}

//---------------------------------------------------------------------------//
void sample_sum( rtt_dsxx::UnitTest &ut, bool const omrpn )
{
    if( rtt_c4::node() == 0 )
        std::cout << "Begin test sample_sum()...\n" << std::endl;
    
    // Generate data and benchmark values:
    int N(10000000);
    std::vector<double> foo(N,0.0);
    std::vector<double> result(N,0.0);
    std::vector<double> bar(N,99.0);

    Timer t1_serial_build;
    t1_serial_build.start();
    
    for(int i=0;i<N;++i)
    {
        foo[i] = 99.00+i;
        bar[i] = 0.99*i;
        result[i] = std::sqrt(foo[i]+bar[i])+1.0;
    }
    t1_serial_build.stop();
    
    Timer t2_serial_accumulate;
    t2_serial_accumulate.start();
    
    double sum = std::accumulate(foo.begin(), foo.end(), 0.0);

    t2_serial_accumulate.stop();
    
    if( node() == 0 )
        std::cout << "benchmark: sum(foo) = " << sum << std::endl;

#ifdef USE_OPENMP
    {
        // More than 1 MPI rank per node --> turn off OMP.
        if( ! omrpn )
            omp_set_num_threads( 1 );

        if( node() == 0 )
            std::cout << "\nNow computing sum using " << omp_get_max_threads()
                      << " OMP threads." << std::endl;
        
        // Generate omp_result 
        std::vector<double> omp_result(N,0.0);
        double omp_sum(0.0);

        Timer t1_omp_build;
        t1_omp_build.start();

#pragma omp parallel for shared(foo,bar)
        for(int i=0;i<N;++i)
        {
            foo[i] = 99.00+i;
            bar[i] = 0.99*i;
            result[i] = std::sqrt(foo[i]+bar[i])+1.0;
        }
        t1_omp_build.stop();
        
        // Accumulate via OMP

        Timer t2_omp_accumulate;
        t2_omp_accumulate.start();
        
#pragma omp parallel for reduction(+: omp_sum)
        for( int i=0; i<N; ++i )
            omp_sum += foo[i];

        t2_omp_accumulate.stop();
        
        // Sanity check
        if( rtt_dsxx::soft_equiv(sum,omp_sum) )
            PASSMSG( "OpenMP sum matches std::accumulate() value!" );
        else
            FAILMSG( "OpenMP sum differs!" );

        if( node() == 0 )
        {
            std::cout.precision(6);
            std::cout.setf(std::ios::fixed,std::ios::floatfield);
            std::cout << "Timers:"
                      << "\n\t             \tSerial Time \tOMP Time"
                      << "\n\tbuild      = \t" << t1_serial_build.wall_clock()
                      << "\t" << t1_omp_build.wall_clock()
                      << "\n\taccumulate = \t" << t2_serial_accumulate.wall_clock()
                      << "\t" << t2_omp_accumulate.wall_clock() << std::endl;
        }
        
        if( omrpn )
        {
            if( t2_omp_accumulate.wall_clock()
                < t2_serial_accumulate.wall_clock() )
                PASSMSG( "OMP accumulate was faster than Serial accumulate.");
            else
                FAILMSG( "OMP accumulate was slower than Serial accumulate.");
        }        
    }
#endif
    return;
}

//---------------------------------------------------------------------------//
// This is a simple demonstration problem for OMP.  Nothing really to check
// for PASS/FAIL. 
int MandelbrotCalculate(complex c, int maxiter)
{
    // iterates z = z*z + c until |z| >= 2 or maxiter is reached, returns the
    // number of iterations

    complex z = c;
    int n = 0;
    for(; n<maxiter; ++n)
    {
        if( std::abs(z) >= 2.0 ) break;
        z = z*z+c;
    }
    return n;
}

void MandelbrotDriver(rtt_dsxx::UnitTest & ut)
{
    if( rtt_c4::node() == 0 )
        std::cout << "\nGenerating Mandelbrot image (OMP threads)...\n"
                  << std::endl;
    
    const int width  = 78;
    const int height = 44;
    const int num_pixels = width*height;
    const complex center(-0.7, 0.0);
    const complex span(   2.7, -(4/3.0)*2.7*height/width );
    const complex begin = center-span/2.0;
    // const complex end   = center+span/2.0;
    const int maxiter = 100000;

    // Use OMP threads
    Timer t;
    t.start();

#pragma omp parallel for ordered schedule(dynamic)
    for( int pix=0; pix<num_pixels; ++pix)
    {
        const int x = pix%width;
        const int y = pix/width;

        complex c = begin + complex( x * span.real() / (width+1.0),
                                     y * span.imag() / (height+1.0) );

        int n = MandelbrotCalculate(c,maxiter);
        if( n == maxiter) n = 0;

#pragma omp ordered
        {
            char c = ' ';
            if( n>0 )
            {
                static const char charset[] = ".,c8M@jawrpogOQEPGJ";
                c = charset[ n% (sizeof(charset)-1) ];
            }
            std::putchar(c);
            if(x+1 == width) std::puts("|");
        }
    }
    t.stop();
    double const gen_time_omp = t.wall_clock();

    // Repeat for serial case
    if( rtt_c4::node() == 0 )
        std::cout << "\nGenerating Mandelbrot image (Serial)...\n"
                  << std::endl;
    
    t.reset();
    t.start();

    for( int pix=0; pix<num_pixels; ++pix)
    {
        const int x = pix%width;
        const int y = pix/width;

        complex c = begin + complex( x * span.real() / (width+1.0),
                                     y * span.imag() / (height+1.0) );

        int n = MandelbrotCalculate(c,maxiter);
        if( n == maxiter) n = 0;

        {
            char c = ' ';
            if( n>0 )
            {
                static const char charset[] = ".,c8M@jawrpogOQEPGJ";
                c = charset[ n% (sizeof(charset)-1) ];
            }
            std::putchar(c);
            if(x+1 == width) std::puts("|");
        }
    }
    t.stop();
    double const gen_time_serial = t.wall_clock();

    std::cout << "\nTime to generate Mandelbrot:"
              << "\n   Normal: " << gen_time_serial << " sec."
              << "\n   OMP   : " << gen_time_omp    << " sec." << std::endl;

    if( gen_time_omp < gen_time_serial )
        PASSMSG( "OMP generation of Mandelbrot image is faster.");
    else
        FAILMSG( "OMP generation of Mandelbrot image is slower.");
    
    return;
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        // One MPI rank per machine node?
        bool omrpn(false);
        
        // Unit tests
        topo_report( ut, omrpn );
        sample_sum( ut, omrpn );
        if( rtt_c4::nodes() == 1 )
            MandelbrotDriver(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstOMP, " 
                  << err.what() << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstOMP, " 
                  << "An unknown exception was thrown." << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of phw.cc
//---------------------------------------------------------------------------//
