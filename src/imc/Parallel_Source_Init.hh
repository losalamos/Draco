//----------------------------------*-C++-*----------------------------------//
// Parallel_Source_Init.hh
// Todd J. Urbatsch
// Mon Aug  3 09:31:56 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __imc_Parallel_Source_Init_hh__
#define __imc_Parallel_Source_Init_hh__

//===========================================================================//
// class Parallel_Source_Init - 
//
// Purpose : Parallel_Source_Init calculates source energies on multiple
//           processors and collects the information at the master node.
//           The master node iterates on number of particles of each source 
//           type, then distributes the processor-dependent source numbers
//           to all processors.
//
//
// revision history:
// -----------------
//  0) original - only for Full Domain Decomposition. (Requires procs_per_cell
//                which could be assembled from cells_per_proc and which
//                divides numbers of particles to each processor and the
//                random number stream number.  Some generality exists in
//                the accumulation of energies.)
//  1) 032699 : readded the T^4 slopes that are provided to the interface by
//              the host code
// 
//===========================================================================//

#include "Opacity.hh"
#include "Mat_State.hh"
#include "Particle.hh"
#include "Particle_Buffer.hh"
#include "Source.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <string>
#include <vector>
#include <iostream>

namespace rtt_imc 
{

template<class MT, class PT = Particle<MT> >
class Parallel_Source_Init 
{
  private:
    // data received from MT_Interface
    std::vector<double> evol_ext;
    std::vector<double> rad_source;
    double rad_s_tend;
    std::vector<std::string> ss_pos;
    std::vector<double> ss_temp;
    std::vector<double> rad_temp;
    std::vector<int> global_cell;
    double delta_t;
    double elapsed_t;
    int npmax;
    int npnom;
    double dnpdt;
    std::string ss_dist;
    int num_global_cells;
    int cycle;
    
    // source initialization data

    // number of particles for this cycle
    int npwant;

    // volume source variables
    typename MT::CCSF_double evol;
    typename MT::CCSF_double evol_net;
    double evoltot;

    // surface source variables
    typename MT::CCSF_double ess;
    typename MT::CCSF_int fss;
    double esstot;

    // radiation energy per cell, total for census energy
    typename MT::CCSF_double ecen;
    double ecentot;

    // number of census particles per cell
    typename MT::CCSF_int ncen;
    int ncentot;
    dsxx::SP<typename Particle_Buffer<PT>::Census> census;

    // number of surface source and volume source particles
    typename MT::CCSF_int nvol;
    typename MT::CCSF_int nss;
    int nvoltot;
    int nsstot;

    // energy loss due to inadequate sampling of evol, ss, and initial census
    double eloss_vol;
    double eloss_ss;
    double eloss_cen;

    // energy weights for census, ss, and vol emission source particles
    typename MT::CCSF_double ew_vol;
    typename MT::CCSF_double ew_ss;
    typename MT::CCSF_double ew_cen;

    // random number streams for each source
    typename MT::CCSF_int volrn;
    typename MT::CCSF_int ssrn;
    typename MT::CCSF_int cenrn;

    // T^4 slope for source
    typename MT::CCVF_double t4_slope;

    // maximum number of cells capable of fitting on a processor
    int capacity;

    // global source vectors (only on master node)
    std::vector<double> global_ecen;
    std::vector<double> global_evol;
    std::vector<double> global_ess;
    std::vector<int> global_ncen;
    std::vector<int> global_nvol;
    std::vector<int> global_nss;
    std::vector<double> global_ew_cen;
    std::vector<double> global_ew_vol;
    std::vector<double> global_ew_ss;
    std::vector<int> global_cenrn;
    std::vector<int> global_volrn;
    std::vector<int> global_ssrn;

    // global source energies and losses
    double global_ecentot;
    double global_evoltot;
    double global_esstot;
    int    global_ncentot;
    int    global_nvoltot;
    int    global_nsstot;
    double global_eloss_cen;
    double global_eloss_vol;
    double global_eloss_ss;

    // cells on each processor, only used on the master processor
    std::vector<std::vector<int> > cells_on_proc;

    // number of source particles, census, source energies, number of volume
    // and surface sources
    void calc_source_energies(const Opacity<MT> &, const Mat_State<MT> &);
    void calc_source_numbers(const Opacity<MT> &);
    void comb_census(const MT &, rtt_rng::Rnd_Control &);

    // initial census service functions
    void calc_evol(const Opacity<MT> &, const Mat_State<MT> &);
    void calc_ess();
    void calc_init_ecen();
    void sum_up_ecen(const MT &);
    void calc_ncen_init();
    void write_initial_census(const MT &, rtt_rng::Rnd_Control &);

    // update the global total number of census particles after combing
    void update_global_ncentot();

    // calculate slope of T_electron^4 for volume emission
    void calc_t4_slope(const MT &, const Mat_State<MT> &);
    
    // communication functions
    void send_source_energies(const MT &);
    void recv_source_energies(const MT &);
    void send_source_numbers(const MT &);
    void recv_source_numbers(const MT &);
    void send_census_numbers(const MT &);
    void recv_census_numbers(const MT &);

    // typedefs
    typedef dsxx::SP<typename Particle_Buffer<PT>::Census> Census_SP;

  public:
    // constructor for master node
    template<class IT> Parallel_Source_Init(dsxx::SP<IT>, dsxx::SP<MT>);

    // initial census function
    Census_SP calc_initial_census(dsxx::SP<MT>, dsxx::SP<Opacity<MT> >,
				  dsxx::SP<Mat_State<MT> >, 
				  dsxx::SP<rtt_rng::Rnd_Control>);

    // source initialyzer function
    dsxx::SP<Source<MT> > initialize(dsxx::SP<MT>, dsxx::SP<Opacity<MT> >,  
				     dsxx::SP<Mat_State<MT> >, 
				     dsxx::SP<rtt_rng::Rnd_Control>, 
				     const Particle_Buffer<PT> &); 

    // accessor functions
    typename MT::CCSF_double get_evol_net() const { return evol_net; }
    inline int get_global_ncentot() const;
    inline int get_global_nsstot() const;
    inline int get_global_nvoltot() const;

    // diagnostic functions
    void print(std::ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT, class PT>
inline std::ostream& operator<<(std::ostream &out, 
				const Parallel_Source_Init<MT,PT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
// get total (global) number of census particles

template<class MT, class PT>
int Parallel_Source_Init<MT,PT>::get_global_ncentot() const
{
    // only call on the Master node
    Require (!node());
    return global_ncentot;
}

//---------------------------------------------------------------------------//
// get total (global) number of surface source particles

template<class MT, class PT>
int Parallel_Source_Init<MT,PT>::get_global_nsstot() const
{
    // only call on the Master node
    Require (!node());
    return global_nsstot;
}

//---------------------------------------------------------------------------//
// get total (global) number of volume emission particles

template<class MT, class PT>
int Parallel_Source_Init<MT,PT>::get_global_nvoltot() const
{
    // only call on the Master node
    Require (!node());
    return global_nvoltot;
}


} // end namespace rtt_imc

#endif                          // __imc_Parallel_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Parallel_Source_Init.hh
//---------------------------------------------------------------------------//
