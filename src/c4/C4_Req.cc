//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Req.cc
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req member definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "C4_Req.hh"
#include "ds++/Assert.hh"
#include <iostream>

namespace rtt_c4 {

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//---------------------------------------------------------------------------//

C4_Req::C4_Req() : p(new C4_ReqRefRep) { ++p->n; }

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 *
 * Attach to an existing message request.
 */
//---------------------------------------------------------------------------//

C4_Req::C4_Req(const C4_Req &req) : p(NULL) {
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

C4_Req::~C4_Req() { free_(); }

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment.
 *
 * Detach from our prior message request, waiting on it if necessary.  Then
 * attach to the new one.
 */
//---------------------------------------------------------------------------//

C4_Req &C4_Req::operator=(const C4_Req &req) {
  free_();

  if (req.inuse())
    p = req.p;
  else
    p = new C4_ReqRefRep;

  ++p->n;

  return *this;
}

//---------------------------------------------------------------------------//
/*!
 * Utility for cleaning up letter in letter/envelope idiom
 */
/* private */
void C4_Req::free_() {
  --p->n;
  if (p->n <= 0)
    delete p;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * Register a new non blocking message request.
 */
//---------------------------------------------------------------------------//

C4_ReqRefRep::C4_ReqRefRep()
    : n(0), assigned(0)
#ifdef C4_MPI
      ,
      s(MPI_Status()), r(MPI_Request())
#endif
{
  // empty
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 *
 * It is important that all existing requests are cleared before the
 * destructor is called.  We used to have a wait() in here; however, this
 * causes exception safety problems.  In any case, it is probably a bad idea
 * to clean up communication by going out of scope.
 */
//---------------------------------------------------------------------------//

C4_ReqRefRep::~C4_ReqRefRep() { /* empty */
}

//---------------------------------------------------------------------------//
//! Wait for an asynchronous message to complete.
//---------------------------------------------------------------------------//

void C4_ReqRefRep::wait() {
  if (assigned) {
#ifdef C4_MPI
    MPI_Wait(&r, &s);
#endif
  }
  clear();
}

//---------------------------------------------------------------------------//
//! Tests for the completion of a non blocking operation.
//---------------------------------------------------------------------------//

bool C4_ReqRefRep::complete() {
#ifdef C4_MPI
  int flag = 0;
  bool indicator = false;
  if (assigned)
    MPI_Test(&r, &flag, &s);
  if (flag != 0) {
    clear();
    Check(r == MPI_REQUEST_NULL);
    indicator = true;
  }
  return indicator;
#endif
#ifdef C4_SCALAR
  throw "Send to self machinery has not been implemented in scalar mode.";
#endif
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of C4_Req.cc
//---------------------------------------------------------------------------//
