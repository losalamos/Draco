//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle.t.hh"
 * \author Thomas M. Evans, Todd J. Urbatsch, Mike Buksas
 * \date   Fri Jan 30 17:04:24 1998
 * \brief  Particle class implementation file.
 *
 * This file is included in the Particle.hh header file so that the Particle
 * base class does not have to be explicitly instantiated.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

//===========================================================================//
// PARTICLE<MT,FT> FUNCTIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// STATIC MEMBER VARIABLES
//---------------------------------------------------------------------------//
/*!
 *
 * \brief Parameter minwt_frac is the fractional energy-weight cutoff between
 * implicit and analog absorption behavior.
 */
template<class MT> const double Particle<MT>::minwt_frac = 0.01;

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Convert a particle event descriptor string into an int.
 */
template<class MT>
int Particle<MT>::convert_string_to_descriptor(std_string desc)
{
    // declare return type
    int return_value;

    // born descriptors
    if (desc == "born")
	return_value = BORN;
    else if (desc == "census_born")
	return_value = CENSUS_BORN;
    else if (desc == "boundary_born")
	return_value = BOUNDARY_BORN;
    else if (desc == "vol_emission")
	return_value = VOL_EMISSION;
    else if (desc == "surface_source")
	return_value = SURFACE_SOURCE;
    else if (desc == "unpacked")
	return_value = UNPACKED;

    // collision event descriptors
    else if (desc == "scatter")
	return_value = SCATTER;
    else if (desc == "low_weight")
	return_value = LOW_WEIGHT;
    else if (desc == "eff_scatter")
	return_value = EFF_SCATTER;
    else if (desc == "thom_scatter")
	return_value = THOM_SCATTER;
 
    // streaming descriptors
    else if (desc == "reflection")
	return_value = REFLECTION;
    else if (desc == "stream")
	return_value = STREAM;
    else if (desc == "escape")
	return_value = ESCAPE;
    else if (desc == "cross_boundary")
	return_value = CROSS_BOUNDARY;

    // time and census descriptors
    else if (desc == "census")
	return_value = CENSUS;

    // death
    else if (desc == "killed")
	return_value = KILLED;

    // last else
    else 
	Insist(0,"Invalid string descriptor");

    // return
    return return_value;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert an int into a particle event descriptor.
 */
template<class MT>
std::string Particle<MT>::convert_descriptor_to_string(int index)
{
    std_string desc;

    switch (index) {

	// born descriptors
    case BORN:           
	desc = "born"; 
	break;
    case CENSUS_BORN:    
	desc = "census_born"; 
	break;
    case BOUNDARY_BORN:  
	desc = "boundary_born"; 
	break;
    case VOL_EMISSION:  
	desc = "vol_emission"; 
	break;
    case SURFACE_SOURCE: 
	desc = "surface_source"; 
	break;
    case UNPACKED:
	desc = "unpacked";
	break;

	// collision event descriptors
    case SCATTER:       
	desc = "scatter";
	break;
    case LOW_WEIGHT: 	 
	desc = "low_weight";
	break;
    case EFF_SCATTER:	
	desc = "eff_scatter";
	break;
    case THOM_SCATTER:
	desc = "thom_scatter";
	break;
	
	// streaming descriptors
    case REFLECTION:	
	desc = "reflection";
	break;
    case STREAM:	 
	desc = "stream";
	break;
    case ESCAPE:	
	desc = "escape";
	break;
    case CROSS_BOUNDARY:
	desc = "cross_boundary";
	break;
	
	// time and census descriptors
    case CENSUS:	 
	desc =  "census";
	break;
	
	// death
    case KILLED:	 
	desc =  "killed";
	break;
	
    default:
	Insist(0,"Unrecognized descriptor number");
	break;
    }  

    return desc;
}

//---------------------------------------------------------------------------//
// DIAGNOSTIC MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out a particle to some stream.
 */
template<class MT>
void Particle<MT>::print(std::ostream &output) const
{
    using std::ios;
    using std::setiosflags;
    using std::endl;
    using std::setw;

    // set precisions
    output.precision(3);
    output << setiosflags(ios::fixed);
    
    output << "*** PARTICLE DATA ***" << endl; 
    output << "---------------------" << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < r.size(); i++)
	output << setw(12) << r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < omega.size(); i++)
	output << setw(12) << omega[i] << " ";
    output << endl;
    
    // cell
    output << setw(20) << setiosflags(ios::right) << "Cell: " << setw(12) 
	   << cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Energy-weight: " 
           << setw(12) << ew << endl;
}

//===========================================================================//
// CLASS PARTICLE<MT>::DIAGNOSTIC
//===========================================================================//

template<class MT>
void Particle<MT>::Diagnostic::header() const 
{ 
    output << "*** PARTICLE HISTORY ***" << std::endl; 
    output << "------------------------" << std::endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print(const Particle<MT> &particle)
    const
{
    using std::ios; 

    // set output precision
    output.precision(3);
    output.setf(ios::scientific, ios::floatfield);

    // print particulars of the particle based on its status
    if (particle.alive == true)
	print_alive(particle);
    else
	print_dead(particle);
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_alive(const Particle<MT> &particle) const 
{
    using std::endl;
    using std::ios;
    using std::setiosflags;
    using std::setw;

    // print active particle (alive = true)
    output << " -- Particle is alive -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << convert_descriptor_to_string(particle.descriptor) 
	   << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;
    
    // cell
    output << setw(20) << setiosflags(ios::right) << "Cell: " << setw(12) 
	   << particle.cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Energy-weight: " 
           << setw(12) << particle.ew << endl;

    // fraction of original weight
    output << setw(20) << setiosflags(ios::right) << "Fraction: " 
           << setw(12) << particle.fraction << endl;

    // time remaining in this time step
    output << setw(20) << setiosflags(ios::right) << "Time_Left: " 
           << setw(12) << particle.time_left << endl;
    
    output << endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_dead(const Particle<MT> &particle) const
{
    using std::endl;
    using std::ios;
    using std::setw;
    using std::setiosflags;

    // print dead particle (alive = false)
    output << " -- Particle is dead -- " << endl;

    // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << convert_descriptor_to_string(particle.descriptor) 
	   << endl;
    
    // coordinates
    output << setw(20) << setiosflags(ios::right) << " Last Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
    // direction
    output << setw(20) << setiosflags(ios::right) << " Last Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;

    // cell
    output << setw(20) << setiosflags(ios::right) << "Last Cell: " 
	   << setw(12) << particle.cell << endl;
    
    // energy-weight, ew
    output << setw(20) << setiosflags(ios::right) << "Last Energy-weight: "
	   << setw(12) << particle.ew << endl;

    // fraction of original weight
    output << setw(20) << setiosflags(ios::right) << "Last Fraction: " 
           << setw(12) << particle.fraction << endl;

    // time remaining in this time step
    output << setw(20) << setiosflags(ios::right) << "Last Time_Left: " 
           << setw(12) << particle.time_left << endl;

    output << endl;
}

//---------------------------------------------------------------------------//

template<class MT>
void Particle<MT>::Diagnostic::print_dist(double d_scat, 
					  double d_bnd, 
					  double d_cen, 
					  int    cell) const
{
    using std::ios;
    using std::setw;
    using std::endl;
    using std::setiosflags;

    // do detailed diagnostic print of particle event distances
    output << setw(20) << setiosflags(ios::right) << "Present cell: "
	   << setw(12) << cell << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-scatter: "
	   << setw(12) << d_scat << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-boundary: "
	   << setw(12) << d_bnd << endl;   
    output << setw(20) << setiosflags(ios::right) << "Dist-census: "
	   << setw(12) << d_cen << endl;
}

//---------------------------------------------------------------------------//
//                              end of Particle.t.hh
//---------------------------------------------------------------------------//
