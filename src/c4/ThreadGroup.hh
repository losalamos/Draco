//----------------------------------*-C++-*----------------------------------//
// ThreadGroup.hh
// Geoffrey M. Furnish
// Mon Jul 13 12:56:20 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_ThreadGroup_hh__
#define __c4_ThreadGroup_hh__

#include <pthread.h>

#include "config.hh"

#include <vector>

C4_NAMESPACE_BEG

class ThreadControl;

//===========================================================================//
// class ThreadGroup - Manage a collection of cooperating threads

// This class is used to manage a set of cooperating threads.  Creation of a
// thread group results in spawning a collection of POSIX threads.  An
// instance the template parameter class is created in each thread, and run.
// This is kind of like the Java Runnable concept of threading.
//===========================================================================//

template<class T>
class ThreadGroup {

    int nthreads;

    ThreadControl *tcb;

    std::vector<T*> v;

    ThreadGroup( const ThreadGroup& );
    ThreadGroup& operator=( const ThreadGroup& );
  public:
    ThreadGroup( int nthreads_, const T& model = T() );
    ~ThreadGroup();
};

C4_NAMESPACE_END

#endif                          // __c4_ThreadGroup_hh__

//---------------------------------------------------------------------------//
//                              end of c4/ThreadGroup.hh
//---------------------------------------------------------------------------//
