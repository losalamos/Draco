/*-----------------------------------*-C-*-----------------------------------*/
/* imc/Source_Builder_Defs.h
 * Thomas M. Evans 
 * Thu Aug  1 11:50:07 2002 */
/*---------------------------------------------------------------------------*/
/* $Id$ */
/*---------------------------------------------------------------------------*/

#include <imc/config.h>

/* Scoping rules required in the Source_Builder inheritance tree
   using the COMPAQ CXX compiler */

#ifdef IMC_COMPAQ

/* Protected variables */
#define npnom            Source_Builder<MT,FT,PT>::npnom
#define npmax            Source_Builder<MT,FT,PT>::npmax
#define dnpdt            Source_Builder<MT,FT,PT>::dnpdt
#define cycle            Source_Builder<MT,FT,PT>::cycle
#define delta_t          Source_Builder<MT,FT,PT>::delta_t
#define ss_dist          Source_Builder<MT,FT,PT>::ss_dist
#define census           Source_Builder<MT,FT,PT>::census
#define topology         Source_Builder<MT,FT,PT>::topology
#define parallel_data_op Source_Builder<MT,FT,PT>::parallel_data_op
#define npwant           Source_Builder<MT,FT,PT>::npwant
#define ecen             Source_Builder<MT,FT,PT>::ecen  
#define ecentot          Source_Builder<MT,FT,PT>::ecentot
#define ew_cen           Source_Builder<MT,FT,PT>::ew_cen
#define evol             Source_Builder<MT,FT,PT>::evol
#define evol_net         Source_Builder<MT,FT,PT>::evol_net
#define mat_vol_src      Source_Builder<MT,FT,PT>::mat_vol_src
#define evoltot          Source_Builder<MT,FT,PT>::evoltot
#define ew_vol           Source_Builder<MT,FT,PT>::ew_vol
#define mat_vol_srctot   Source_Builder<MT,FT,PT>::mat_vol_srctot
#define ess              Source_Builder<MT,FT,PT>::ess
#define ss_face_in_cell  Source_Builder<MT,FT,PT>::ss_face_in_cell
#define esstot           Source_Builder<MT,FT,PT>::esstot
#define ew_ss            Source_Builder<MT,FT,PT>::ew_ss
#define volrn            Source_Builder<MT,FT,PT>::volrn
#define ssrn             Source_Builder<MT,FT,PT>::ssrn
#define freq_samp_data   Source_Builder<MT,FT,PT>::freq_samp_data

/* Protected functions */
#define calc_source_energies   Source_Builder<MT,FT,PT>::calc_source_energies
#define calc_initial_ecen      Source_Builder<MT,FT,PT>::calc_initial_ecen
#define calc_num_src_particles Source_Builder<MT,FT,PT>::calc_num_src_particles 
#define write_initial_census   Source_Builder<MT,FT,PT>::write_initial_census 
#define comb_census            Source_Builder<MT,FT,PT>::comb_census 
#define reset_ew_in_census     Source_Builder<MT,FT,PT>::reset_ew_in_census

#endif

/*---------------------------------------------------------------------------*/
/*                              end of imc/Source_Builder_Defs.h */
/*---------------------------------------------------------------------------*/
