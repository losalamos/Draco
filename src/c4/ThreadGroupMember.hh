//----------------------------------*-C++-*----------------------------------//
// ThreadGroupMember.hh
// Geoffrey M. Furnish
// Tue Jul 14 11:21:45 1998
//---------------------------------------------------------------------------//
// @> Declare class representing thread family participation.
//---------------------------------------------------------------------------//

#ifndef __c4_ThreadGroupMember_hh__
#define __c4_ThreadGroupMember_hh__

#include "config.hh"
#include "ThreadControl.hh"

#include "ds++/Assert.hh"

#include <pthread.h>

C4_NAMESPACE_BEG

// Scalar reduction applicative template classes, ala expression template
// kernels.  Note that the init() methods don't have to do anything in the
// current formulation of the scalar_reduce routine, since it is just setting 
// the address of the reduction point to be the same as the address of one of 
// the participating scalars, so it's effectively already initialized.

template<class T>
class ScalarSum {
  public:
    static void init ( T& result, const T& x ) {}
    static void apply( T& result, const T& x ) { result += x; }
};

template<class T>
class ScalarProd {
  public:
    static void init ( T& result, const T& x ) {}
    static void apply( T& result, const T& x ) { result *= x; }
};

template<class T>
class ScalarMin {
  public:
    static void init ( T& result, const T& x ) {}
    static void apply( T& result, const T& x ) 
    {
        if (x < result) result = x;
    }
};

template<class T>
class ScalarMax {
  public:
    static void init ( T& result, const T& x ) {}
    static void apply( T& result, const T& x ) 
    {
        if (x > result) result = x;
    }
};

//===========================================================================//
// class ThreadGroupMember - Represents participation in a family of threads

// This class serves as a base class for objects which will be used to
// implement threaded operations in a cooperative context.  A class may
// publicly derive from ThreadGroupMember to indicate that it "is-a" member
// of a thread group.  
//
// Several services are provided.  First, data is provided to indicate the
// thread's relationship to other threads in its group.  This is by analogy
// with "C4::NodeInfo".  Additionally, we provide a variety of thread
// reduction support routines to facilitate cooperative multithreading
// activities. 
//===========================================================================//

class ThreadGroupMember {

  protected:
    pthread_t tid;

    int thread;
    int threads;

    ThreadControl *tcb;

    template<class T, class ScalarOp> 
    void scalar_reduce( T& x ) const;

  public:
    ThreadGroupMember();
    virtual ~ThreadGroupMember() {}

    virtual void c4_run_thread() =0;

    template<class X> friend class ThreadGroup;

    void gsync() const;

    template<class T> 
    void gsum ( T& x ) const { scalar_reduce< T, ScalarSum<T> >( x ); }
    template<class T> 
    void gprod( T& x ) const { scalar_reduce< T, ScalarProd<T> >( x ); }
    template<class T> 
    void gmin ( T& x ) const { scalar_reduce< T, ScalarMin<T> >( x ); }
    template<class T> 
    void gmax ( T& x ) const { scalar_reduce< T, ScalarMax<T> >( x ); }

    template<class T>
    void gsum( T *px, int n ) const;

    int lock()   const { return pthread_mutex_lock( &tcb->mutex ); }
    int unlock() const { return pthread_mutex_unlock( &tcb->mutex ); }

};

C4_NAMESPACE_END

#include "ThreadGroupMember.t.cc"

#endif                          // __c4_ThreadGroupMember_hh__

//---------------------------------------------------------------------------//
//                              end of c4/ThreadGroupMember.hh
//---------------------------------------------------------------------------//
