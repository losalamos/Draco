/*-----------------------------------*-C-*-----------------------------------*/
/* mc/Unset_Particle_Buffer_Defs.h
 * Thomas M. Evans 
 * Thu Aug  1 10:44:05 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <mc/config.h>

/* Clear scoping rules required in the Particle_Buffer inheritance tree using
   the COMPAQ CXX compiler */

#ifdef MC_COMPAQ

/* Protected variables */
#undef packed_particles
#undef int_data        
#undef comm_int        
#undef comm_char       

/* Public functions */
#undef is_empty                  
#undef get_size_packed_particle  
#undef get_maximum_num_particles 
#undef empty_buffer              

#endif

/*---------------------------------------------------------------------------*/
/*                 end of mc/Unset_Particle_Buffer_Defs.h */
/*---------------------------------------------------------------------------*/
