//----------------------------------*-C++-*----------------------------------//
// C4_Req.hh
// Geoffrey Furnish
// Wed Feb  1 12:58:32 1995
//---------------------------------------------------------------------------//
// @> A class for managing non blocking message requests.
//---------------------------------------------------------------------------//

#ifndef __c4_C4_Req_hh__
#define __c4_C4_Req_hh__

#include "c4/config.hh"

// autodoc: noprint C4_ReqRefRep

C4_NAMESPACE_BEG

struct C4_ReqRefRep {
    int n;

    C4_ReqRefRep() { n = 1; }
};

//===========================================================================//
// class C4_Req - Handle for non blocking message requests

// This class provides an encapsulator for the message id's (NX) or requests
// (MPI) which are produced by non blocking calls.  This class automatically
// waits for the message to complete when the containing object goes out of
// scope, thus plugging one of the easiest types of programming errors with
// non blocking messaging.  Reference counting is used so that these may be
// passed by value without accidentally triggering a program stall.
//===========================================================================//

class C4_Req {

    C4_ReqRefRep *p;

    int assigned;

  public:			// So C4 api functions can see them, but
				// don't anybody else touch them!
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

  public:
    C4_Req();
    C4_Req( const C4_Req& req );
    ~C4_Req();
    C4_Req& operator=( const C4_Req& req );

    void set()      { assigned = 1; }
    void clear()    { assigned = 0; }

    void wait();
    void free();

    bool complete();

    int inuse() const { return assigned; }
};

C4_NAMESPACE_END

#endif                          // __c4_C4_Req_hh__

//---------------------------------------------------------------------------//
//                              end of c4/C4_Req.hh
//---------------------------------------------------------------------------//
