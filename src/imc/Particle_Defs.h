/*-----------------------------------*-C-*-----------------------------------*/
/* imc/Particle_Defs.h
 * Thomas M. Evans 
 * Fri Aug  2 09:27:20 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <imc/config.h>

/* Scoping rules required in the Particle inheritance tree using the COMPAQ
   CXX compiler */

#ifdef IMC_COMPAQ

/* Protected variables */
#define minwt_frac     Particle<MT>::minwt_frac
#define ew             Particle<MT>::ew
#define r              Particle<MT>::r
#define omega          Particle<MT>::omega
#define cell           Particle<MT>::cell
#define time_left      Particle<MT>::time_left
#define fraction       Particle<MT>::fraction
#define alive          Particle<MT>::alive
#define descriptor     Particle<MT>::descriptor
#define random         Particle<MT>::random

/* Protected functions */
#define stream                Particle<MT>::stream
#define stream_analog_capture Particle<MT>::stream_analog_capture
#define surface               Particle<MT>::surface
#define use_analog_absorption Particle<MT>::use_analog_absorption
#define census_event          Particle<MT>::census_event
#define boundary_event        Particle<MT>::boundary_event
#define scatter               Particle<MT>::scatter

/* Protected enumeration */
#define BORN           Particle<MT>::BORN
#define CENSUS_BORN    Particle<MT>::CENSUS_BORN
#define BOUNDARY_BORN  Particle<MT>::BOUNDARY_BORN
#define VOL_EMISSION   Particle<MT>::VOL_EMISSION
#define SURFACE_SOURCE Particle<MT>::SURFACE_SOURCE
#define UNPACKED       Particle<MT>::UNPACKED
#define SCATTER        Particle<MT>::SCATTER
#define LOW_WEIGHT     Particle<MT>::LOW_WEIGHT
#define EFF_SCATTER    Particle<MT>::EFF_SCATTER
#define THOM_SCATTER   Particle<MT>::THOM_SCATTER 
#define COLLISION      Particle<MT>::COLLISION
#define CUTOFF         Particle<MT>::CUTOFF
#define REFLECTION     Particle<MT>::REFLECTION
#define STREAM         Particle<MT>::STREAM
#define ESCAPE         Particle<MT>::ESCAPE
#define CROSS_BOUNDARY Particle<MT>::CROSS_BOUNDARY
#define ABSORPTION     Particle<MT>::ABSORPTION
#define BOUNDARY       Particle<MT>::BOUNDARY
#define CENSUS         Particle<MT>::CENSUS
#define KILLED         Particle<MT>::KILLED

/* Diagnostic sub-class */
#define output Particle<MT>::Diagnostic::output

/* Public functions */
#define status Particle<MT>::status

#endif

/*---------------------------------------------------------------------------*/
/*                              end of imc/Particle_Defs.h */
/*---------------------------------------------------------------------------*/
