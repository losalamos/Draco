/*-----------------------------------*-C-*-----------------------------------*/
/* LAMG/c2f.h
 * Randy M. Roberts 
 * Tue Jan 25 14:30:23 2000 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#ifndef __LAMG_c2f_h__
#define __LAMG_c2f_h__

#include <LAMG/config.h>

#if CCALLSFSUB_CASE == CCALLSFSUB_LOWER

#if CCALLSFSUB_METHOD == CCALLSFSUB_SUFFIX
#define CCALLFSUB(NUP,NLOW) NLOW##_
#elif CCALLSFSUB_METHOD == CCALLSFSUB_PREFIX
#define CCALLFSUB(NUP,NLOW) _##NLOW
#elif CCALLSFSUB_METHOD == CCALLSFSUB_NONE
#define CCALLFSUB(NUP,NLOW) NLOW
#else
#error Invalid CCALLSFSUB_METHOD
#endif

#elif CCALLSFSUB_CASE == CCALLSFSUB_UPPER

#if CCALLSFSUB_METHOD == CCALLSFSUB_SUFFIX
#define CCALLFSUB(NUP,NLOW) NUP##_
#elif CCALLSFSUB_METHOD == CCALLSFSUB_PREFIX
#define CCALLFSUB(NUP,NLOW) _##NUP
#elif CCALLSFSUB_METHOD == CCALLSFSUB_NONE
#define CCALLFSUB(NUP,NLOW) NUP
#else
#error Invalid CCALLSFSUB_METHOD
#endif

#else

#error Invalid CCALLSFSUB_CASE

#endif

#define CCALLSFSUB0(NUP,NLOW) \
  CCALLFSUB(NUP,NLOW)()
#define CCALLSFSUB1(NUP,NLOW,T1,A1) \
  CCALLFSUB(NUP,NLOW)(A1)
#define CCALLSFSUB2(NUP,NLOW,T1,A1,T2,A2) \
  CCALLFSUB(NUP,NLOW)(A1,A2)
#define CCALLSFSUB3(NUP,NLOW,T1,A1,T2,A2,T3,A3) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3)
#define CCALLSFSUB4(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4)
#define CCALLSFSUB5(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4,T5,A5) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4,A5)
#define CCALLSFSUB6(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4,T5,A5,T6,A6) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4,A5,A6)
#define CCALLSFSUB7(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4,T5,A5,T6,A6,T7,A7) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4,A5,A6,A7)
#define CCALLSFSUB8(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4,T5,A5,T6,A6,T7,A7,T8,A8) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4,A5,A6,A7,A8)
#define CCALLSFSUB9(NUP,NLOW,T1,A1,T2,A2,T3,A3,T4,A4,T5,A5, \
                             T6,A6,T7,A7,T8,A8,T9,A9) \
  CCALLFSUB(NUP,NLOW)(A1,A2,A3,A4,A5,A6,A7,A8,A9)

#endif                          /* __LAMG_c2f_h__ */

/*---------------------------------------------------------------------------*/
/*                              end of LAMG/c2f.h */
/*---------------------------------------------------------------------------*/
