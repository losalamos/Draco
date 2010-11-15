//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstglobal_containers.cc
 * \author Kent Budge
 * \date   Mon Mar 24 09:41:04 2008
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <set>

#include "ds++/Assert.hh"
#include "../Release.hh"
#include "c4/ParallelUnitTest.hh"
#include "c4/C4_Functions.hh"
#include "../global_containers.i.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstglobal_containers(UnitTest &ut)
{
    unsigned const pid = rtt_c4::node();
    unsigned const number_of_processors = rtt_c4::nodes();

    {
        set<unsigned> local_set;
        local_set.insert(pid);
        local_set.insert(number_of_processors+pid);
        
        global_merge(local_set);
        
        if (local_set.size()==2*number_of_processors)
        {
            ut.passes("Correct number of global elements");
        }
        else
        {
            ut.failure("NOT correct number of global elements");
        }
        
        for (unsigned p=0; p<number_of_processors; ++p)
        {
            if (local_set.count(p)!=1 ||
                local_set.count(number_of_processors+p)!=1)
            {
                ut.failure("WRONG element in set");
            }
        }
    }

    {
        map<unsigned, double> local_map;
        local_map[pid] = pid;
        local_map[number_of_processors+pid] = 2*pid;
        
        global_merge(local_map);
        
        if (local_map.size()==2*number_of_processors)
        {
            ut.passes("Correct number of global elements");
        }
        else
        {
            ut.failure("NOT correct number of global elements");
        }
        
        for (unsigned p=0; p<number_of_processors; ++p)
        {
            if (local_map.count(p)!=1 ||
                local_map.count(number_of_processors+p)!=1)
            {
                ut.failure("WRONG element in map");
            }
            if (local_map[p]!=p ||
                local_map[number_of_processors+p]!=2*p)
            {
                ut.failure("WRONG element value in map");
            }
        }
    }

    {
        map<unsigned, bool> local_map;
        local_map[pid] = false;
        local_map[number_of_processors+pid] = true;
        
        global_merge(local_map);
        
        if (local_map.size()==2*number_of_processors)
        {
            ut.passes("Correct number of global elements");
        }
        else
        {
            ut.failure("NOT correct number of global elements");
        }
        
        for (unsigned p=0; p<number_of_processors; ++p)
        {
            if (local_map.count(p)!=1 ||
                local_map.count(number_of_processors+p)!=1)
            {
                ut.failure("WRONG element in map");
            }
            if (local_map[p]!=false ||
                local_map[number_of_processors+p]!=true)
            {
                ut.failure("WRONG element value in map");
            }
        }
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstglobal_containers(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstglobal_containers, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstglobal_containers, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstglobal_containers.cc
//---------------------------------------------------------------------------//
