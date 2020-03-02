//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/C4_Req.cc
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req member definitions.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "C4_Req.hh"

namespace rtt_c4 {

//----------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//----------------------------------------------------------------------------//
C4_Req::C4_Req() : p(new C4_ReqRefRep) { ++p->n; }

//----------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 *
 * Attach to an existing message request.
 */
//----------------------------------------------------------------------------//
C4_Req::C4_Req(const C4_Req &req) : p(NULL) {
  if (req.inuse())
    p = req.p;
  else
    p = new C4_ReqRefRep;
  ++p->n;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * If we've been left holding the bag, make sure the message has completed.
 * This should plug a wide class of potential programming errors.
 */
//----------------------------------------------------------------------------//
C4_Req::~C4_Req() { free_(); }

//----------------------------------------------------------------------------//
/*!
 * \brief Assignment.
 *
 * Detach from our prior message request, waiting on it if necessary.  Then
 * attach to the new one.
 */
//----------------------------------------------------------------------------//
C4_Req &C4_Req::operator=(const C4_Req &req) {
  free_();

  if (req.inuse())
    p = req.p;
  else
    p = new C4_ReqRefRep;

  ++p->n;

  return *this;
}

//----------------------------------------------------------------------------//
/*!
 * Utility for cleaning up letter in letter/envelope idiom
 */
//----------------------------------------------------------------------------//
/* private */
void C4_Req::free_() {
  --p->n;
  if (p->n <= 0)
    delete p;
}

void C4_ReqRefRep::free() {
#ifdef C4_MPI
  if (assigned) {
    MPI_Cancel(&r);
    MPI_Request_free(&r);
  }
#endif
  clear();
}

//----------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//----------------------------------------------------------------------------//
C4_ReqRefRep::C4_ReqRefRep()
    : n(0), assigned(false)
#ifdef C4_MPI
      ,
      r(MPI_Request())
#endif
{
  // empty
}

//----------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * It is important that all existing requests are cleared before the destructor
 * is called.  We used to have a wait() in here; however, this causes exception
 * safety problems.  In any case, it is probably a bad idea to clean up
 * communication by going out of scope.
 */
//----------------------------------------------------------------------------//
C4_ReqRefRep::~C4_ReqRefRep() { /* empty */
}

//----------------------------------------------------------------------------//
/*!
 * \brief Wait for an asynchronous message to complete.
 * \param status Status object.
 *
 * This function is non-const because it updates the underlying request data
 * member.
 */
// ---------------------------------------------------------------------------//
#ifdef C4_MPI

void C4_ReqRefRep::wait(C4_Status *status) {
  if (assigned) {
    MPI_Status *s = MPI_STATUS_IGNORE;
    if (status) {
      s = status->get_status_obj();
      Check(s);
    }
    MPI_Wait(&r, s);
  }
  clear();
}

#elif defined(C4_SCALAR)

void C4_ReqRefRep::wait(C4_Status * /*status*/) { clear(); }

#endif

//----------------------------------------------------------------------------//
/*!
 * \brief Tests for the completion of a non blocking operation.
 * \param status Status object.
 *
 * This function is non-const because it updates the underlying request data
 * member.
 */
//----------------------------------------------------------------------------//
#ifdef C4_MPI

bool C4_ReqRefRep::complete(C4_Status *status) {
  int flag = 0;
  bool indicator = false;
  if (assigned) {
    MPI_Status *s = MPI_STATUS_IGNORE;
    if (status) {
      s = status->get_status_obj();
      Check(s);
    }
    MPI_Test(&r, &flag, s);
  }
  if (flag != 0) {
    clear();
    Check(r == MPI_REQUEST_NULL);
    indicator = true;
  }
  return indicator;
  // throw "C4_Req::complete() has not been implemented for this communicator";
}

#elif defined(C4_SCALAR)

bool C4_ReqRefRep::complete(C4_Status * /*status*/) {
  return true;
  //  throw "C4_Req::complete() has not been implemented in scalar mode.";
}

#endif

} // end namespace rtt_c4

//----------------------------------------------------------------------------//
// end of C4_Req.cc
//----------------------------------------------------------------------------//
