//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Mat_State.t.hh
 * \author Thomas M. Evans
 * \date   Mon Mar  9 16:06:28 1998
 * \brief  Mat_State class member template function definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Mat_State.hh"
#include <iomanip>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
  
 * \brief Mat_State constructor.

 * \param density_ cell-centered field of material densities
 * \param temp_ cell-centered field of material (electron) temperatures
 * \param dedt_ cell-centered field of heat capacities
 * \param spec_heat_ cell-centered field of specific heat capacities

 */
template<class MT>
Mat_State<MT>::Mat_State(const ccsf_double &density_, 
			 const ccsf_double &temp_,
			 const ccsf_double &dedt_,
			 const ccsf_double &spec_heat_)
    : density(density_),
      temperature(temp_), 
      dedt(dedt_), 
      spec_heat(spec_heat_)
{
    Ensure (density.get_Mesh() == temperature.get_Mesh());
    Ensure (density.get_Mesh() == dedt.get_Mesh());
    Ensure (density.get_Mesh() == spec_heat.get_Mesh());
}

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// diagnostic for printing

template<class MT>
void Mat_State<MT>::print(std::ostream &output, int cell) const
{
  // print out material state of cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;
    
    output.precision(4);
    
    output << setw(8) << cell << setw(15) << setiosflags(ios::scientific)
	   << density(cell) << setw(15) << temperature(cell) << setw(15) 
	   << dedt(cell) << setw(15) << spec_heat(cell) << endl;
} 

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

template<class MT>
std::ostream& operator<<(std::ostream &output, const Mat_State<MT> &object)
{
  // print out opacities for all cells
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;
    
    output << setw(8)  << setiosflags(ios::right) << "Cell" 
	   << setw(15) << "Density" << setw(15) << "Temp" 
	   << setw(15) << "dEdT" << setw(15) << "Cv" << endl;
    output << "--------------------------------------"
	   << "------------------------------" << endl;

    for (int i = 1; i <= object.num_cells(); i++)
	object.print(output, i);
    return output;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Mat_State.t.hh
//---------------------------------------------------------------------------//
