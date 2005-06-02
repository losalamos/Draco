//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Req.cc
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req member definitions.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "ds++/Assert.hh"
#include "C4_Req.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//---------------------------------------------------------------------------//

C4_Req::C4_Req()
    : p(new C4_ReqRefRep)
{
    ++p->n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 *
 * Attach to an existing message request.
 */
//---------------------------------------------------------------------------//

C4_Req::C4_Req( const C4_Req& req )
{
    if (req.inuse())
        p = req.p;
    else
        p = new C4_ReqRefRep;
    ++p->n;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * If we've been left holding the bag, make sure the message has completed.
 * This should plug a wide class of potential programming errors.
 */
//---------------------------------------------------------------------------//

C4_Req::~C4_Req()
{
    --p->n;
    if (p->n <= 0)
	delete p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment.
 *
 * Detach from our prior message request, waiting on it if necessary.  Then
 * attach to the new one.
 */
//---------------------------------------------------------------------------//

C4_Req& C4_Req::operator=( const C4_Req& req )
{
    --p->n;
    if (p->n <= 0)
	delete p;

    if (req.inuse())
        p = req.p;
    else
        p = new C4_ReqRefRep;

    ++p->n;

    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//---------------------------------------------------------------------------//

C4_ReqRefRep::C4_ReqRefRep()
    : assigned(0), n(0)
{
    // empty
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * If we've been left holding the bag, make sure the message has completed.
 * This should plug a wide class of potential programming errors.
 */
//---------------------------------------------------------------------------//

C4_ReqRefRep::~C4_ReqRefRep()
{
    wait();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait for an asynchronous message to complete.
 */
//---------------------------------------------------------------------------//

void C4_ReqRefRep::wait()
{
    if (assigned) 
    {
#ifdef C4_MPI
	MPI_Wait(&r, &s);
#endif
    }
    clear();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Free request handle for a posted asynchronous receive.
 *
 * Note: once freed the handle must be reactivated to test for completeness or
 * to wait on it
 */
//---------------------------------------------------------------------------//

void C4_ReqRefRep::free()
{
#ifdef C4_MPI
    if (assigned)
	MPI_Cancel( &r );
#endif
    clear();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tests for the completion of a non blocking operation.
 */
//---------------------------------------------------------------------------//

bool C4_ReqRefRep::complete()
{
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

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of items returned on the last complete operation.
 */
//---------------------------------------------------------------------------//

unsigned C4_ReqRefRep::count()
{
#ifdef C4_MPI
    int count;
    MPI_Get_count( &s, MPI_CHAR, &count );
    return count;
#else
    return 0;
#endif
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                 end of C4_Req.cc
//---------------------------------------------------------------------------//
