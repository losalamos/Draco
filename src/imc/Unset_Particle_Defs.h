/*-----------------------------------*-C-*-----------------------------------*/
/* imc/Particle_Defs.h
 * Thomas M. Evans 
 * Fri Aug  2 09:27:20 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <imc/config.h>

/* Clear scoping rules required in the Particle inheritance tree using the
   COMPAQ CXX compiler */

#ifdef IMC_COMPAQ

/* Protected variables */
#undef minwt_fraction 
#undef ew             
#undef r              
#undef omega          
#undef cell           
#undef time_left      
#undef fraction       
#undef alive          
#undef descriptor     
#undef random         

/* Protected functions */
#undef stream                
#undef stream_analog_capture 
#undef surface               
#undef use_analog_absorption 
#undef census_event          
#undef boundary_event        
#undef scatter               

/* Protected enumeration */
#undef BORN           
#undef CENSUS_BORN    
#undef BOUNDARY_BORN  
#undef VOL_EMISSION   
#undef SURFACE_SOURCE 
#undef UNPACKED       
#undef SCATTER        
#undef LOW_WEIGHT     
#undef EFF_SCATTER    
#undef THOM_SCATTER   
#undef COLLISION      
#undef CUTOFF         
#undef REFLECTION     
#undef STREAM         
#undef ESCAPE         
#undef CROSS_BOUNDARY 
#undef ABSORPTION     
#undef BOUNDARY       
#undef CENSUS         
#undef KILLED 

/* Diagnostic sub-class */
#undef output

/* Public functions */
#undef status     

#endif

/*---------------------------------------------------------------------------*/
/*                              end of imc/Particle_Defs.h */
/*---------------------------------------------------------------------------*/
