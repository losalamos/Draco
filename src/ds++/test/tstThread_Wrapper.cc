//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstThread_Wrapper.cc
 * \author Tim Kelley
 * \date   Aug. 14, 2017
 * \note   Copyright (C) 2006-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

#include "ds++/Thread_Wrapper.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <chrono>
#include <functional>   // std::function
#include <iostream>
#include <sstream>
#include <string>

#define FAIL_IF_NOT(c) if(!(c)) ITFAILS
#define FAIL_IF(c) if((c)) ITFAILS

using rtt_dsxx::Thread_Wrapper;
using rtt_dsxx::UnitTest;

/* Three thread actions: */
void quick_action(){
  std::chrono::milliseconds so_long(1);
  std::this_thread::sleep_for(so_long);
  return;
}

void slow_action1(){
  std::chrono::seconds so_long(1);
  std::this_thread::sleep_for(so_long);
  return;
}

void slow_action2(std::stringstream &s){
  std::chrono::seconds so_long(1);
  std::this_thread::sleep_for(so_long);
  s << "slow_action: done\n";
  return;
}

/**\brief Exercise default ctor */
void testInstantiation0(UnitTest &/*ut*/) {
  Thread_Wrapper tw;
  return;
} // testInstantiation1

/**\brief Exercise two arg ctor with the join action. */
void testInstantiation2(UnitTest &ut) {
  {
    Thread_Wrapper tw(std::thread(quick_action),Thread_Wrapper::action::join);
  }
  // Construting a thread and move it into a Thread_Wrapper
  // NB: cannot copy a std::thread.
  {
    std::thread t(slow_action1);
    auto tid1 = t.get_id();
    Thread_Wrapper tw(std::move(t),Thread_Wrapper::action::join);
    auto tid2 = tw.get().get_id();
    FAIL_IF_NOT(tid1 == tid2);
  }
  return;
} // testInstantiation2

/**\brief Exercise two arg ctor with the detach action.
 *
 * This test is a bit contrived: the calling thread will detach from
 * the worker thread. The worker thread will continue to run, i.e. sleep then
 * write a message to the stream. Meanwhile, the caller will print a message,
 * then sleep. So barring some sort of 2000 msec catastrophe, the writes shd
 * be ordered. If that fails, it's not clear it says anything about
 * Thread_Wrapper; it may just be that something odd happened in the
 * execution. So if this test starts failing, it might need to be retired
 * or modified.
 */
void testDetach(UnitTest &ut){
  std::stringstream s;
  {
    std::thread t(slow_action2,std::ref(s));
    auto tid = t.get_id();
    Thread_Wrapper tw(std::move(t),Thread_Wrapper::action::detach);
    auto tid2 = tw.get().get_id();
    FAIL_IF_NOT(tid == tid2);
  }
  s << "host thread: done\n";
  std::chrono::seconds so_long(1);
  std::this_thread::sleep_for(so_long);
  // It would be odd if the detached OS thread finished faster than its
  // std::thread / Thread_Wrapper handle went out of scope.
  FAIL_IF_NOT(s.str() == "host thread: done\nslow_action: done\n");
  return;
} //testDetach

using t_func = std::function<void(UnitTest &)>;

void run_a_test(UnitTest & u, t_func f, std::string const & msg){
  f(u);
  if (u.numFails == 0){
    u.passes(msg);
  }
  return;
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    run_a_test(ut,testInstantiation0,"Thread_Wrapper::ctor() ok.");
    run_a_test(ut,testInstantiation2,
      "Thread_Wrapper::ctor(std::thread && t, action a) ok.");
    run_a_test(ut,testDetach,
      "Thread_Wrapper::ctor(std::thread && t, detach) ok.");
  } // try--catches in the epilog:
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of ds++/test/tstThread_Wrapper.cc
//---------------------------------------------------------------------------//

