//----------------------------------*-C++-*----------------------------------//
// C4_Req.cc
// Geoffrey Furnish
// Wed Feb  1 12:58:33 1995
//---------------------------------------------------------------------------//
// @> A class for managing non blocking message requests.
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "C4_Req.hh"
#include "ds++/Assert.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Constructor.  Register a new non blocking message request.
//---------------------------------------------------------------------------//

C4_Req::C4_Req()
    : assigned(0)
{
    p = new C4_ReqRefRep;
}

//---------------------------------------------------------------------------//
// Copy constructor.  Attach to an existing message request.
//---------------------------------------------------------------------------//

C4_Req::C4_Req( const C4_Req& req )
{
    assigned = req.assigned;

    p = req.p;
    p->n++;

    PGN( mid = req.mid );
    MPI( r = req.r );
#ifdef C4_SHMEM
    mid = req.mid;
    thread = req.thread;
    type = req.type;
#endif
}

//---------------------------------------------------------------------------//
// Destructor.  If we've been left holding the bag, make sure the message has
// completed.  This should plug a wide class of potential programming errors.
//---------------------------------------------------------------------------//

C4_Req::~C4_Req()
{
    p->n--;
    if (!p->n) {
	delete p;
	wait();
    }
}

//---------------------------------------------------------------------------//
// Assignment.  Detach from our prior message request, waiting on it if
// necessary.  Then attach to the new one.
//---------------------------------------------------------------------------//

C4_Req& C4_Req::operator=( const C4_Req& req )
{
    p->n--;
    if (!p->n) {
	delete p;
	wait();
    }

    assigned = req.assigned;

    p = req.p;
    p->n++;

    PGN( mid = req.mid );
    MPI( r = req.r );
#ifdef C4_SHMEM
    mid = req.mid;
    thread = req.thread;
    type = req.type;
#endif

    return *this;
}

//---------------------------------------------------------------------------//
// Wait for an asynchronous message to complete.
//---------------------------------------------------------------------------//

void C4_Req::wait()
{
    if (assigned) {
	PGN( msgwait( mid ) );
	MPI( MPI_Wait( &r, &s ) );
#ifdef C4_SHMEM
// 	C4_Wait( mid, type );
	C4_Wait( thread );
#endif
    }
    clear();
}

//---------------------------------------------------------------------------//
// Free request handle for a posted asynchronous receive, note: once freed
// the handle must be reactivated to test for completeness or to wait on it
//---------------------------------------------------------------------------//

void C4_Req::free()
{
#ifdef C4_SHMEM
    throw "incomplete";
#endif
#ifdef C4_MPI
    if (assigned)
	MPI_Cancel( &r );
#endif
    clear();
}

//---------------------------------------------------------------------------//
// Tests for the completion of a non blocking operation.
//---------------------------------------------------------------------------//

bool C4_Req::complete()
{
#ifdef C4_SHMEM
    throw "incomplete";
#endif
#ifdef C4_MPI
    int flag       = 0;
    bool indicator = false;
    if (assigned)
        MPI_Test( &r, &flag, &s );
    if (flag != 0)
    {
        clear();
	Check ( r == MPI_REQUEST_NULL);
	indicator = true;
    }
    return indicator;
#endif
#ifdef C4_SCALAR
    throw "Send to self machinery has not been implemented in scalar mode.";
#endif
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of C4_Req.cc
//---------------------------------------------------------------------------//
