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
#include "Global.hh"
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

 */
template<class MT, class PT>
Source<MT, PT>::Source(ccsf_int &vol_rnnum_, 
		       ccsf_int &nvol_,
		       ccsf_double &ew_vol_,
		       ccsf_int &ss_rnnum_, 
		       ccsf_int &nss_, 
		       ccsf_int &fss_,
		       ccsf_double &ew_ss_,
		       SP_Census census_,
		       std::string ssd, 
		       int nvoltot_, 
		       int nsstot_,
		       SP_Rnd_Control rcon_, 
		       SP_Mat_State mat_state,
		       SP_Mesh_Op operations,
		       SP_Topology top) 
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
      topology(top)
{
    Require (topology);

    // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh()      == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh()      == nss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_vol.get_Mesh());
    Check (nvol.get_Mesh()      == fss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_ss.get_Mesh());

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
template<class MT, class PT>
rtt_dsxx::SP<PT> Source<MT, PT>::get_Source_Particle(double delta_t)
{
    using rtt_dsxx::SP;

    bool sampled = false;

    // instantiate particle to return
    SP<PT> source_particle;

    // do all surface source particles, one-by-one
    while (!sampled && nssdone < nsstot)
    {
	if (nsrcdone_cell < nss(current_cell))
	{
	    int rn_str_id =
		rtt_mc::global::mod_with_2e9(ss_rnnum(current_cell) +
					     nsrcdone_cell);  
	    rcon->set_num(rn_str_id);
	    source_particle = get_ss(delta_t);
	    sampled = true;
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
		Check (nssdone == nsstot);
		current_cell = 1;
	    }
	}
    }

    // do all volume emission particles, one-by-one
    while (!sampled && nvoldone < nvoltot)
    {
	if (nsrcdone_cell < nvol(current_cell))
	{
	    int rn_str_id =
		rtt_mc::global::mod_with_2e9(vol_rnnum(current_cell) + 
					     nsrcdone_cell);
	    rcon->set_num(rn_str_id);
	    source_particle = get_evol(delta_t);
	    sampled = true;
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
		Check (nvoldone == nvoltot);
		current_cell = 1;
	    }
	}
    }

    // do all census particles, one-by-one
    while (!sampled && ncendone < ncentot)
    {
	source_particle = get_census(delta_t);
	sampled = true;
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
template<class MT, class PT>
rtt_dsxx::SP<PT> Source<MT, PT>::get_ss(double delta_t)
{
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
	double costheta = sqrt(rand.ran());
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
    SP<PT> ss_particle(new PT(r, omega, ew, cell, rand, fraction, 
			      time_left, "surface_source"));

    // return the ss particle;
    return ss_particle;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a volume source particle.
 */
template<class MT, class PT>
rtt_dsxx::SP<PT> Source<MT, PT>::get_evol(double delta_t)
{
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
    SP<PT> vol_particle(new PT(r, omega, ew, cell, rand, fraction, 
			       time_left, "vol_emission")); 

    return vol_particle;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a census particle.
 */
template<class MT, class PT>
rtt_dsxx::SP<PT> Source<MT, PT>::get_census(double delta_t)
{
    using rtt_dsxx::SP;

    Require (census->size() > 0);

    // get the census particle from the Census buffer
    SP<PT> census_particle = census->top();
    census_particle->set_time_left(delta_t);
    census_particle->set_descriptor("from_census");

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
// DIAGNOSTICS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out the source to an output.

 * \param out ostream output

 */
template<class MT, class PT>
void Source<MT, PT>::print(std::ostream &out) const
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
