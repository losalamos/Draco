//----------------------------------*-C++-*----------------------------------//
// C4_Req.hh
// Geoffrey Furnish
// Wed Feb  1 12:58:32 1995
//---------------------------------------------------------------------------//
// @> A class for managing non blocking message requests.
//---------------------------------------------------------------------------//

#ifndef __c4_C4_Req_hh__
#define __c4_C4_Req_hh__

// C4 package configure
#include <c4/config.h>

#ifdef C4_MPI
#include "c4_mpi.h"
#endif

namespace rtt_c4
{

//===========================================================================//
// class C4_ReqRefRep - Handle for non blocking message requests

// This class provides an encapsulator for the message id's (NX) or requests
// (MPI) which are produced by non blocking calls.  This class automatically
// waits for the message to complete when the containing object goes out of
// scope, thus plugging one of the easiest types of programming errors with
// non blocking messaging.  Reference counting is used so that these may be
// passed by value without accidentally triggering a program stall.
//===========================================================================//

class C4_ReqRefRep {

    friend class C4_Req;
    
    int n;
    int assigned;

#ifdef C4_NX
    long mid;
#endif
#ifdef C4_MPI
    MPI_Status  s;
    MPI_Request r;
#endif
#ifdef C4_SHMEM
    pthread_t thread;
    int mid;
    int type;
#endif

  private:

    // Disallowed methods

    C4_ReqRefRep( const C4_ReqRefRep& rep );
    C4_ReqRefRep& operator=( const C4_ReqRefRep& rep );

    // Private default ctor and dtor for access from C4_Req only.
    
    C4_ReqRefRep();
    ~C4_ReqRefRep();

  public:
    
    void wait();
    void free();

    bool complete();

    int inuse() const { return assigned; }

  private:
    
    void set()      { assigned = 1; }
    void clear()    { assigned = 0; }
};

//===========================================================================//
// class C4_Req - Envelope for non blocking message requests

// This class provides an encapsulator for the message id's (NX) or requests
// (MPI) which are produced by non blocking calls.  This class automatically
// waits for the message to complete when the containing object goes out of
// scope, thus plugging one of the easiest types of programming errors with
// non blocking messaging.  Reference counting is used so that these may be
// passed by value without accidentally triggering a program stall.
//===========================================================================//

class C4_Req {
    
    C4_ReqRefRep *p;

  public:
    
    C4_Req();
    C4_Req( const C4_Req& req );
    ~C4_Req();
    C4_Req& operator=( const C4_Req& req );

    void wait()     { p->wait(); }
    void free()     { p->free(); }

    bool complete() { return p->complete(); }

    int inuse() const { return p->inuse(); }

  private:

    void set()      { p->set(); }
    void clear()    { p->clear(); }

    // Private access to the C4_ReqRefRep internals.

#ifdef C4_MPI
    MPI_Status  &s() { return p->s; }
    MPI_Request &r() { return p->r; }
#endif

    // FRIENDSHIP
    
    // Specific friend C4 functions that may need to manipulate the
    // C4_ReqRefRep internals.

#ifdef C4_MPI
    template<class T>
    friend C4_Req send_async(const T *buf, int nels, int dest, int tag);
    template<class T>
    friend C4_Req receive_async(T *buf, int nels, int source, int tag);
    template<class T>
    friend void send_async(C4_Req& r, const T *buf, int nels, int dest,
			   int tag);
    template<class T>
    friend void receive_async(C4_Req& r, T *buf, int nels, int source, 
			      int tag);
#endif
};

} // end namespace rtt_c4

#endif                          // __c4_C4_Req_hh__

//---------------------------------------------------------------------------//
//                              end of c4/C4_Req.hh
//---------------------------------------------------------------------------//
