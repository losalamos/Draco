/*-----------------------------------*-C-*-----------------------------------*/
/* mc/Particle_Buffer_Defs.h
 * Thomas M. Evans 
 * Thu Aug  1 10:44:05 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <mc/config.h>

/* Scoping rules required in the Particle_Buffer inheritance tree
   using the COMPAQ CXX compiler */

#ifdef MC_COMPAQ

/* Protected variables */
#define packed_particles Particle_Buffer<PT>::packed_particles
#define int_data         Particle_Buffer<PT>::int_data
#define comm_int         Particle_Buffer<PT>::comm_int
#define comm_char        Particle_Buffer<PT>::comm_char

/* Public functions */
#define is_empty                  Particle_Buffer<PT>::is_empty
#define get_size_packed_particle  Particle_Buffer<PT>::get_size_packed_particle
#define get_maximum_num_particles Particle_Buffer<PT>::get_maximum_num_particles
#define empty_buffer              Particle_Buffer<PT>::empty_buffer

#endif

/*---------------------------------------------------------------------------*/
/*                              end of mc/Particle_Buffer_Defs.h */
/*---------------------------------------------------------------------------*/
