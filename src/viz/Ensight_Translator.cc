//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Translator.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator implementation file (non-templated code).
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Ensight_Translator.hh"

namespace rtt_viz
{

using std::endl;
using std::setw;
using std::ios;
using std::ofstream;
using std::string;

//---------------------------------------------------------------------------//
// WRITE ENSIGHT CASE FILE
//---------------------------------------------------------------------------//
/*!
 * \brief Write out case file.
 */
void Ensight_Translator::ensight_case(const double)
{
    // create the case file name (directory already created)
    const char *filename = case_filename.c_str();
    
    // open the case file
    ofstream caseout(filename);

    // write the format header
    caseout << "FORMAT" << endl;
    caseout << "type: ensight" << endl << endl;

    // write the geometry file block
    caseout << "GEOMETRY" << endl;
    caseout << "model: 1   " << "./geo/data.****" << endl << endl;

    // write the variable block header
    caseout << "VARIABLE" << endl;

    // write the pointer to the node variables
    for (int i = 0; i < ens_vdata_names.size(); i++)
	caseout << " scalar per node:    1  " << setw(19)
		<< ens_vdata_names[i] << setw(4) << " "
		<< "./" << ens_vdata_names[i] << "/data.****" << endl;

    // write the pointer to the cell variables
    for (int i = 0; i < ens_cdata_names.size(); i++)
	caseout << " scalar per element: 1  " << setw(19) 
		<< ens_cdata_names[i] << setw(4) << " "
		<< "./" << ens_cdata_names[i] << "/data.****" << endl;

    caseout << endl;
    // write out the time block
    caseout << "TIME" << endl;
    caseout << "time set:              " << setw(4) << "   1" << endl;
    caseout << "number of steps:       " << setw(4) << dump_times.size() 
	    << endl;
    caseout << "filename start number: " << setw(4) << "   1" << endl;
    caseout << "filename increment:    " << setw(4) << "   1" << endl;
    caseout << "time values:           " << endl;
    
    // write out times
    caseout.precision(5);
    caseout.setf(ios::scientific, ios::floatfield);
    for (int i = 0; i < dump_times.size(); i++)
	caseout << setw(12) << dump_times[i] << endl;
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                              end of Ensight_Translator.cc
//---------------------------------------------------------------------------//
