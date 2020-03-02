//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   plot2D/plot2D_grace.h
 * \author Thomas M. Evans
 * \brief  Solve a system of equations by Gaussian elimination.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef __plot2D_plot2D_grace_h__
#define __plot2D_plot2D_grace_h__

#include "ds++/Assert.hh"
#include "plot2D/config.h"

#ifdef GRACE_H

// ... then Grace is supported on this platform.

/* include grace headers */
#include GRACE_H

bool rtt_plot2D::Plot2D::is_supported() { return true; }

#else

// ... then Grace is not supported on this platform.

bool rtt_plot2D::Plot2D::is_supported() { return false; }

// Mirror the grace functions to avoid link errors.

int GraceOpenVA(char * /*exe*/, int /*bs*/, ...) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceOpen(int /*bs*/) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceIsOpen(void) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceClose(void) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceClosePipe(void) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceFlush(void) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GracePrintf(const char *, ...) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

int GraceCommand(const char *) {
  Insist(0, "Serious Plot2D error.");
  return 1;
}

#endif // GRACE_H

#endif /* __plot2D_plot2D_grace_h__ */

/*---------------------------------------------------------------------------*/
/* end of plot2D/plot2D_grace.h */
/*---------------------------------------------------------------------------*/
