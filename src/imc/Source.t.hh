//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source.t.hh
 * \author Thomas M. Evans, Todd J. Urbatsch
 * \date   Thu May 14 08:45:49 1998
 * \brief  Source class template definitions file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Source.hh"
#include "Gray_Particle.hh"
#include "Multigroup_Particle.hh"
#include "mc/Sampler.hh"
#include "ds++/Assert.hh"
#include <iomanip>
#include <cmath>
#include <vector>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for source class.

 * \param vol_rnnum_ cell-centered, scalar field of beginning random number
 * id for volume emission source per cell

 * \param nvol_ cell-centered, scalar field of number of volume emission
 * sources per cell

 * \param ew_vol_ cell-centered, scalar field of energy weights of volume
 * emission sources per cell

 * \param ss_rnnum_ cell-centered, scalar field of beginning random number id
 * for surface source per cell

 * \param nss_ cell-centered, scalar field of number of surface sources per
 * cell

 * \param fss_ cell-centered, scalar field of surface source face position
 * per cell

 * \param ew_ss_ cell-centered, scalar field of energy weights of surface
 * sources per cell

 * \param census_ census bank for this source

 * \param ssd surface source distribution string

 * \param nvoltot_ number of volume emission sources in this source

 * \param nsstot_ number of surface sources in this source

 * \param rcon_ rtt_rng::Rnd_Control object

 * \param mat_state material state object

 * \param operations mesh operations object (controls T^4 slope sampling)

 * \param top topology

 * \param opacity_ smart pointer to the opacity

 * \param data Frequency_Sampling_Data object that hold the probabilities for
 * sampling a straight Planckian when FT=Multigroup_Frequency (ie. the
 * calculation is multigroup)

 */
template<class MT, class FT, class PT>
Source<MT,FT,PT>::Source(ccsf_int                             &vol_rnnum_, 
			 ccsf_int                             &nvol_,
			 ccsf_double                          &ew_vol_,
			 ccsf_int                             &ss_rnnum_, 
			 ccsf_int                             &nss_, 
			 ccsf_int                             &fss_,
			 ccsf_double                          &ew_ss_,
			 SP_Census                             census_,
			 std::string                           ssd, 
			 int                                   nvoltot_, 
			 int                                   nsstot_,
			 SP_Rnd_Control                        rcon_, 
			 SP_Mat_State                          mat_state,
			 SP_Mesh_Op                            operations,
			 SP_Topology                           top,
			 SP_Opacity                            opacity_,
			 const Frequency_Sampling_Data<MT,FT> &data) 
    : vol_rnnum(vol_rnnum_),
      nvol(nvol_), 
      ew_vol(ew_vol_), 
      ss_rnnum(ss_rnnum_), 
      nss(nss_), 
      fss(fss_), 
      ew_ss(ew_ss_),
      ss_dist(ssd), 
      census(census_), 
      nvoltot(nvoltot_), 
      nsstot(nsstot_),
      ncentot(census->size()), 
      rcon(rcon_),
      material(mat_state), 
      mesh_op(operations),
      topology(top),
      opacity(opacity_),
      freq_samp_data(data)
{
    Require (topology);
    Require (material);
    Require (opacity);

    // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh()      == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh()      == nss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_vol.get_Mesh());
    Check (nvol.get_Mesh()      == fss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_ss.get_Mesh());
    Check (opacity->num_cells() == material->num_cells());

    // nsrcdone_cell is the running number of source particles completed for
    // a particular source type in a particular cell.
    nsrcdone_cell = 0;

    // Begin with first cell
    current_cell = 1;

    // running totals of completed source particles, by type
    nssdone  = 0;
    nvoldone = 0;
    ncendone = 0;
}

//---------------------------------------------------------------------------//
// MAIN INTERFACE
//---------------------------------------------------------------------------//
/*!

 * \brief Return a rtt_dsxx::SP to the next emitted source particle (PT)

 * This member creates a particle and returns it to the client.  It keeps
 * track of the number of particles of each type created.  It throws a DBC
 * exception if a particle is asked for after all particles have been
 * created.

 * The particles are created with "active" status.

 * \param delta_t problem timestep

 * \return source particle

 */
template<class MT, class FT, class PT>
rtt_dsxx::SP<PT> Source<MT,FT,PT>::get_Source_Particle(double delta_t)
{
    using rtt_dsxx::SP;

    // sampled flag
    bool sampled = false;

    // instantiate particle to return
    SP<PT> source_particle;

    // do all surface source particles, one-by-one
    while (!sampled && nssdone < nsstot)
    {
	if (nsrcdone_cell < nss(current_cell))
	{
	    int rn_str_id = INTEGER_MODULO_1E9(ss_rnnum(current_cell) +
					       nsrcdone_cell);  
	    rcon->set_num(rn_str_id);
	    source_particle = get_ss(delta_t);
	    sampled         = true;
	    nsrcdone_cell++;
	    nssdone++;
	    if (nssdone == nsstot)
	    {
		current_cell  = 1;
		nsrcdone_cell = 0;
	    }
	}
	else
	{
	    current_cell++;
	    nsrcdone_cell = 0;
	    if (current_cell > nss.get_Mesh().num_cells())
	    {
		Insist (nssdone == nsstot,
			"Missed some surface source particles.");
		current_cell = 1;
	    }
	}
    }

    // do all volume emission particles, one-by-one
    while (!sampled && nvoldone < nvoltot)
    {
	if (nsrcdone_cell < nvol(current_cell))
	{
	    int rn_str_id = INTEGER_MODULO_1E9(vol_rnnum(current_cell) + 
					       nsrcdone_cell);
	    rcon->set_num(rn_str_id);
	    source_particle = get_evol(delta_t);
	    sampled         = true;
	    nsrcdone_cell++;
	    nvoldone++;
	    if (nvoldone == nvoltot)
	    {
		current_cell  = 1;
		nsrcdone_cell = 0;
	    }
	}
	else
	{
	    current_cell++;
	    nsrcdone_cell = 0;
	    if (current_cell > nvol.get_Mesh().num_cells())
	    {
		Insist (nvoldone == nvoltot,
			"Missed some volume source particles.");
		current_cell = 1;
	    }
	}
    }

    // do all census particles, one-by-one
    while (!sampled && ncendone < ncentot)
    {
	source_particle = get_census(delta_t);
	sampled         = true;
	ncendone++;
    }

    Insist (sampled, 
	    "Tried to create a source particle after source is empty!");
    Ensure (source_particle->status());
    
    // return the particle
    return source_particle;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Create a surface source particle.
 */
template<class MT, class FT, class PT>
rtt_dsxx::SP<PT> Source<MT,FT,PT>::get_ss(double delta_t)
{
    using rtt_imc::global::Type_Switch;
    using rtt_mc::global::pi;
    using rtt_rng::Sprng;
    using rtt_dsxx::SP;
    using std::vector;

    // face on which surface source resides
    int face = fss(current_cell);

    // get the random number object
    Sprng rand = rcon->get_rn();

    // sample location
    vector<double> r = nss.get_Mesh().sample_pos_on_face
	(current_cell, face, rand);

    // find inward normal, sample direction, and add

    // normal distributed surface source
    vector<double> omega = nss.get_Mesh().get_normal_in(current_cell, face); 
    
    // cosine distribution of surface source about normal if requested, else
    // retain the inward normal to give a "normal" distribution
    if (ss_dist == "cosine")
    {
	double costheta = std::sqrt(rand.ran());
	double phi      = 2.0 * pi * rand.ran();
	nss.get_Mesh().get_Coord().calc_omega(costheta, phi, omega);
    }
    else 
	Insist(ss_dist == "normal", 
	       "Surface source angle distrib is neither cosine nor normal!"); 

    // complete description of surface source particle    
    double ew        = ew_ss(current_cell);
    int cell         = current_cell;
    double fraction  = 1.0;
    double time_left = rand.ran() * delta_t;

    // instantiate particle to return
    SP<PT> ss_particle = make_ss_particle<Type_Switch<PT>::Type>(
	Type_Switch<FT>(), Type_Switch<PT>(), r, omega, ew, cell, rand, 
	fraction, time_left);

    // return the ss particle;
    return ss_particle;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a volume source particle.
 */
template<class MT, class FT, class PT>
rtt_dsxx::SP<PT> Source<MT,FT,PT>::get_evol(double delta_t)
{
    using rtt_imc::global::Type_Switch;
    using rtt_rng::Sprng;
    using rtt_dsxx::SP;
    using std::vector;

    // get the random number object
    Sprng rand = rcon->get_rn();

    // sample location using tilt
    double T         = material->get_T(current_cell);
    vector<double> r = mesh_op->sample_pos_tilt(current_cell, T, rand);

    // sample particle direction
    vector<double> omega = nvol.get_Mesh().get_Coord().
	sample_dir("isotropic", rand); 
 
    double ew        = ew_vol(current_cell);
    int cell         = current_cell;
    double fraction  = 1.0;
    double time_left = rand.ran() * delta_t;

    // instantiate particle to return
    SP<PT> vol_particle = make_vol_particle<Type_Switch<PT>::Type>(
	Type_Switch<FT>(), Type_Switch<PT>(), r, omega, ew, cell, rand, 
	fraction, time_left);

    return vol_particle;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a census particle.
 */
template<class MT, class FT, class PT>
rtt_dsxx::SP<PT> Source<MT,FT,PT>::get_census(double delta_t)
{
    using rtt_dsxx::SP;

    Require (census->size() > 0);

    // get the census particle from the Census buffer
    SP<PT> census_particle = census->top();
    census_particle->set_time_left(delta_t);
    census_particle->set_descriptor(PT::CENSUS_BORN);

    // remove the census particle from the bank
    census->pop();

    // make sure the particle is active
    census_particle->reset_status();

    // give the census particle a local cell index
    int global_cell = census_particle->get_cell();
    int local_cell  = topology->local_cell(global_cell);
    Ensure (local_cell);
    
    census_particle->set_cell(local_cell);

    // return the particle
    return census_particle;
}

//---------------------------------------------------------------------------//
// PRIVATE PARTIAL SPECIALIZATIONS ON FREQUENCY AND PARTICLE TYPE
//---------------------------------------------------------------------------//
// Specialization of make_ss_particle for Gray_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
Source<MT,FT,PT>::SP_Gray_PT Source<MT,FT,PT>::make_ss_particle(
    Switch_Gray,
    Switch_Gray_PT,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   cell,
    const rtt_rng::Sprng &random,
    const double          fraction,
    const double          time_left)
{
    using rtt_dsxx::SP;

    SP_Gray_PT particle(new Gray_PT(r, omega, ew, cell, random, fraction, 
				    time_left, PT::SURFACE_SOURCE));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// Specialization of make_ss_particle for Multigroup_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
Source<MT,FT,PT>::SP_MG_PT Source<MT,FT,PT>::make_ss_particle(
    Switch_MG,
    Switch_MG_PT,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   cell,
    const rtt_rng::Sprng &random,
    const double          fraction,
    const double          time_left)
{
    using rtt_mc::sampler::sample_planckian_frequency;
    using rtt_dsxx::SP;

    // sample group from Planckian
    int group   = 0;
    int counter = 0;
    double freq = 0.0;
    double tss  = freq_samp_data.ss_temperature(cell);

    while (group == 0)
    {
	freq   = sample_planckian_frequency(random, tss);
	group  = opacity->get_Frequency()->find_group_given_a_freq(freq);

	// advance counter
	counter++;

	Insist (counter < 100, "Unable to sample group.");
    }

    Check (group >  0);
    Check (group <= opacity->get_Frequency()->get_num_groups());

    SP_MG_PT particle(new MG_PT(r, omega, ew, cell, random, group,
				fraction, time_left, PT::SURFACE_SOURCE));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// Specialization of make_vol_particle for Gray_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
Source<MT,FT,PT>::SP_Gray_PT Source<MT,FT,PT>::make_vol_particle(
    Switch_Gray,
    Switch_Gray_PT,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   cell,
    const rtt_rng::Sprng &random,
    const double          fraction,
    const double          time_left)
{
    using rtt_dsxx::SP;

    SP_Gray_PT particle(new Gray_PT(r, omega, ew, cell, random, fraction,
				    time_left, PT::VOL_EMISSION));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// Specialization of make_vol_particle for Multigroup_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
Source<MT,FT,PT>::SP_MG_PT Source<MT,FT,PT>::make_vol_particle(
    Switch_MG,
    Switch_MG_PT,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   cell,
    const rtt_rng::Sprng &random,
    const double          fraction,
    const double          time_left)
{
    using rtt_mc::sampler::sample_planckian_frequency;
    using rtt_mc::sampler::sample_bin_from_discrete_cdf;
    using rtt_dsxx::SP;

    // sample to determine if we should sample from the Planckian or sample
    // from sigma * Planckian (opacity-weighted Planckian)
    double ran = random.ran();
    int group  = 0;
    
    // sample from straight Planckian
    if (ran < freq_samp_data.prob_of_straight_Planck_emission(cell))
    {
	// sample group from Planckian
	int counter = 0;
	double freq = 0.0;

	while (group == 0)
	{
	    freq   = sample_planckian_frequency(random, material->get_T(cell));
	    group  = opacity->get_Frequency()->find_group_given_a_freq(freq);

	    // advance counter
	    counter++;

	    Insist (counter < 100, "Unable to sample group.");
	}
    }

    // sample from opacity weighted Planckian cdf
    else
    {
	// get the group index (add 1 since return value is in [0,G-1])
	group = 1 + sample_bin_from_discrete_cdf(
	    random, opacity->get_emission_group_cdf(cell));

	Check (group > 0);
	Check (group <= opacity->get_Frequency()->get_num_groups());
    }

    // make the particle
    SP_MG_PT particle(new MG_PT(r, omega, ew, cell, random, group,
				fraction, time_left, PT::VOL_EMISSION));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// DIAGNOSTICS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out the source to an output.

 * \param out ostream output

 */
template<class MT, class FT, class PT>
void Source<MT,FT,PT>::print(std::ostream &out) const
{
    using std::ios;
    using std::setiosflags;
    using std::setw;
    using std::endl;

    out << endl;
    out << ">>> SOURCE DATA <<<" << endl;
    out << "===================" << endl;

    // numbers of each type of particle
    out << setw(25) << setiosflags(ios::right) << "Census Particles: "
	<< setw(10) << ncentot << endl;
    out << setw(25) << setiosflags(ios::right) << "Volume Particles: "
	<< setw(10) << nvoltot << endl;
    out << setw(25) << setiosflags(ios::right) << "Surface Particles: "
	<< setw(10) << nsstot << endl;

    // lets look at the number of particles in each cell
    out << endl << " ** Sources **" << endl;
    out << setw(10) << " " << setw(15) << setiosflags(ios::right)
	<< "Volume" << setw(5) << " " << setw(15) << setiosflags(ios::right)
	<< "Surface" << endl;
    out << setw(10) << setiosflags(ios::right) << "Cell"
	<< setw(10) << setiosflags(ios::right) << "Number"
	<< setw(10) << setiosflags(ios::right) << "Stream ID"
	<< setw(10) << setiosflags(ios::right) << "Number"
	<< setw(10) << setiosflags(ios::right) << "Stream ID" << endl;
    for (int i = 1; i <= nvol.get_Mesh().num_cells(); i++)
	out << setw(10) << i << setw(10) << nvol(i) 
	    << setw(10) << vol_rnnum(i)  << setw(10) << nss(i)
	    << setw(10) << ss_rnnum(i)   << endl;

    // lets look at the energy weight in each cell
    out << endl << " ** Source Energy ** " << endl;
    out << setw(10) << setiosflags(ios::right) << "Cell" 
	<< setw(15) << setiosflags(ios::right) << "Volume ew" 
	<< setw(15) << setiosflags(ios::right) << "Surface ew" << endl;
    out.precision(3);
    out.setf(ios::scientific, ios::floatfield);
    for (int i = 1; i <= nvol.get_Mesh().num_cells(); i++)
	out << setw(10) << i << setw(15) << ew_vol(i) << setw(15)
	    << ew_ss(i) << endl;
}
	
} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Source.t.hh
//---------------------------------------------------------------------------//
