//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/Procmon.cc
 * \author Kelly Thompson
 * \date   Monday, Apr 22, 2013, 10:10 am
 * \brief  Procmon class for printing runtime system diagnostics (free memory
 *         per node, etc).
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Procmon.hh 5523 2010-11-30 01:12:12Z kellyt $
//---------------------------------------------------------------------------//

#include "Procmon.hh"
#include "c4/C4_Functions.hh"
#include "ds++/path.hh"
#include <fstream>
#include <sstream>

#ifdef USE_PROCMON
// see 'man getrusage'
#include <sys/time.h>
#include <sys/resource.h>

// see 'man getpid'
#include <sys/types.h>
#include <unistd.h>
#endif

namespace rtt_diagnostics
{

//---------------------------------------------------------------------------------------//
std::vector<std::string> tokenize(
    std::string const & source,
    char        const * delimiter_list,
    bool                keepEmpty)
{
    std::vector<std::string> results;

    size_t prev = 0;
    size_t next = 0;

    while ((next = source.find_first_of(delimiter_list, prev)) != std::string::npos)
    {
        if (keepEmpty || (next - prev != 0))
            results.push_back(source.substr(prev, next - prev));
        prev = next + 1;
    }

    if (prev < source.size())
        results.push_back(source.substr(prev));

    return results;
}

//---------------------------------------------------------------------------------------//
/*
 * Notes:  In Cassio's procmon routines, all absolute values are given in
 * bytes instead of kB.  Cassio's procmon report has the following fields:
 *
 * Estimate: mpi_sum_node( 2 * rss_max ) / max_mem_node
 * RSS     : mpi_sum_node( vmrss       ) / max_mem_node
 * Virt    : mpi_sum_node( vmpeak      ) / max_mem_node
 * RSS_MAX : mpi_sum_node( rss_max     ) / max_mem_node
 *
 * where
 * - rss_max is read from getrusage.ru_maxrss;
 * - vmrss and vmpeak are ready from /proc/<pid>/status
 * - max_mem_node is read from /proc/meminfo as MemTotal
 */
void procmon_resource_print( std::string const & identifier,
                             std::ostream & msg )
{
    // ----------------------------------------
    // Resources of interest.
    // ----------------------------------------

    // int proc_ppid(-1);
    std::string proc_name;

    double MemTotal(-1.0);
    //double MemFree(-1.0);

    // int proc_uid (-1);
    int proc_vmpeak(-1);  // The peak size of the virtual memory allocated to the process
    //int proc_vmsize(-1);  // The size of the virtual memory allocated to the process
    int proc_vmrss(-1);   // The amount of memory mapped in RAM ( instead of swapped out )
    // int proc_vmdata(-1);  // The size of the Data segment
    // int proc_vmstk(-1);   // The stack size.

    // ----------------------------------------
    // Find the PID
    // ----------------------------------------

    long long proc_pid(0);
#ifdef USE_PROCMON
    proc_pid = getpid();
#endif
    std::string proc_pid_string;
    {
        std::ostringstream buf;
        buf << proc_pid;
        proc_pid_string = buf.str();
    }

    // ----------------------------------------
    // Examine /proc/meminfo for total memory and free memory.
    // ----------------------------------------
    
    std::string file_meminfo( "/proc/meminfo" );
    Insist( rtt_dsxx::fileExists( file_meminfo ),
        "Could not open /proc/meminfo!  Is this Linux?" );

    std::ifstream fs;
    fs.open( file_meminfo.c_str(), std::fstream::in );
    Check( fs.good() );
    
    while( !fs.eof() )
    {
        std::string line;
        std::getline( fs, line );
        
        // tokenize the string
        std::vector<std::string> tokens = tokenize( line, " \t" );
        if( tokens.size() > 1 ) // protect against empty line
        {
            if( tokens[0] == std::string("MemTotal:") )  MemTotal = atof(tokens[1].c_str()); // kB
            // if( tokens[0] == std::string("MemFree:") )   MemFree  = atof(tokens[1].c_str()); // kB
        }
    }
    fs.close();
    
    Check( MemTotal > 0 );
    // Check( MemFree > 0 );

    // ----------------------------------------
    // Examine /proc/<PID>/status for memory used
    // ----------------------------------------

    std::string file_status_pid( std::string("/proc/") + proc_pid_string
                                 + std::string("/status") );
#ifdef USE_PROCMON
    Insist( rtt_dsxx::fileExists( file_status_pid ),
            std::string("Could not open") + file_status_pid
            + std::string("!  Is this Linux?") );

    fs.open( file_status_pid.c_str(), std::fstream::in );
    Check( fs.good() );

    while( !fs.eof() )
    {
        std::string line;
        std::getline( fs, line );
        // std::cout << line << std::endl;

        // tokenize the string
        std::vector<std::string> tokens = tokenize( line, " \t" );

        if( tokens.size() > 1 ) // protect against empty line.
        {
            if( tokens[0] == std::string("Name:") )   proc_name = tokens[1];
            // if( tokens[0] == std::string("PPid:") )   proc_ppid = atoi(tokens[1].c_str());
            // if( tokens[0] == std::string("Uid:") )    proc_uid  = atoi(tokens[1].c_str());
            if( tokens[0] == std::string("VmPeak:") ) proc_vmpeak = atoi(tokens[1].c_str());
            // if( tokens[0] == std::string("VmSize:") ) proc_vmsize = atoi(tokens[1].c_str());
            if( tokens[0] == std::string("VmRSS:") )  proc_vmrss  = atoi(tokens[1].c_str());
            // if( tokens[0] == std::string("VmData:") ) proc_vmdata = atoi(tokens[1].c_str());
            // if( tokens[0] == std::string("VmStk:") )  proc_vmstk  = atoi(tokens[1].c_str());
        }

    }
    fs.close();
#endif
    
    // ----------------------------------------
    // Use rusage to obtain rss_max
    // ----------------------------------------

    double proc_vmrssmax(-1.0);
#ifdef USE_PROCMON
    struct rusage ruse_getrusage; /* getrusage */
    
    getrusage( RUSAGE_SELF, &ruse_getrusage );
    proc_vmrssmax = ruse_getrusage.ru_maxrss; /* looks to be in kb */
#endif
    
    // ----------------------------------------
    // Print a report
    // ----------------------------------------

    msg // << "\nMemory status:"
        << "[" << rtt_c4::node() << "] " 
        << proc_name << " (pid: " << proc_pid << ")::" << identifier
//              << "\nUid  : " << proc_uid
        << "\tVmPeak : " << proc_vmpeak << " kB (" << proc_vmpeak/MemTotal*100.0 << "%)"
//        << "\tVmSize : " << proc_vmsize << " kB "
        << "\tVmRss : " << proc_vmrss  << " kB (" << proc_vmrss/MemTotal*100.0 << "%)"
        << "\tVmRss_max : " << proc_vmrssmax  << " kB (" << proc_vmrssmax/MemTotal*100.0 << "%)"
        // << "\nVmData : " << proc_vmdata << " kB"
        // << "\nVmStk  : " << proc_vmstk  << " kB"
        // << "\nMemFree: " << MemFree     << " kB"
        << std::endl;
       
//    std::cout << "\nMemTotal: " << MemTotal << " kB"
//              << "\nMemFree : " << MemFree  << " kB"
//               << std::endl;
    
        
    return;
}

} // end namespace rtt_diagnostics

//---------------------------------------------------------------------------//
// end of diagnostics/Procmon.cc
//---------------------------------------------------------------------------//
