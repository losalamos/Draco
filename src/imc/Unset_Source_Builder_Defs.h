/*-----------------------------------*-C-*-----------------------------------*/
/* imc/Unset_Source_Builder_Defs.h
 * Thomas M. Evans 
 * Thu Aug  1 11:50:07 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <imc/config.h>

/* Clear scoping rules required in the Source_Builder inheritance tree
   using the COMPAQ CXX compiler */

#ifdef IMC_COMPAQ

/* Protected variables */
#undef npnom            
#undef npmax            
#undef dnpdt            
#undef cycle            
#undef delta_t          
#undef ss_dist          
#undef census           
#undef topology         
#undef parallel_data_op 
#undef npwant           
#undef ecen             
#undef ecentot          
#undef ew_cen           
#undef evol             
#undef evol_net         
#undef mat_vol_src      
#undef evoltot          
#undef ew_vol           
#undef mat_vol_srctot   
#undef ess              
#undef ss_face_in_cell  
#undef esstot           
#undef ew_ss            
#undef volrn            
#undef ssrn             
#undef freq_samp_data   

/* Protected functions */
#undef calc_source_energies   
#undef calc_initial_ecen      
#undef calc_num_src_particles 
#undef write_initial_census   
#undef comb_census            
#undef reset_ew_in_census     

#endif

/*---------------------------------------------------------------------------*/
/*                              end of imc/Unset_Source_Builder_Defs.h */
/*---------------------------------------------------------------------------*/
