//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Global_Timer.cc
 * \author Kent G. Budge
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Define methods of class Global_Timer, a POSIX standard timer.
 * \note   Copyright (C) 2002-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Timer.hh 7075 2013-04-01 22:48:15Z kellyt $
//---------------------------------------------------------------------------//

#include <iostream>

#include "Global_Timer.hh"

namespace rtt_c4
{
using namespace std;

bool Global_Timer::active_ = false;
     
//---------------------------------------------------------------------------------------//
Global_Timer::Global_Timer(char const *name)
    :
    name_(name)
{
    Require(name != NULL);

    Ensure(name == this->name());
}

//---------------------------------------------------------------------------------------//
/* static */
void Global_Timer::set_activity(bool const active)
{
    if (rtt_c4::node()==0)
    {
        active_ = active;
        
        cout << "***** Global timers are now ";
        if (active)
            cout << "ACTIVE";
        else
            cout << "INACTIVE";
        
        cout << endl;
    }
}

//---------------------------------------------------------------------------------------//
Global_Timer::~Global_Timer()
{
    if (active_)
    {
        cout << endl;
        cout << "Timing report for timer " << name_ << ':' << endl;
        print(cout);
        cout << endl;
    }
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                              end of c4/Global_Timer.cc
//---------------------------------------------------------------------------//
