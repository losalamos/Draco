//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/Hex_Format.cc
 * \author John McGhee
 * \date   Tue Mar  7 08:38:04 2000
 * \brief  Implements a CIC-19 Hex Mesh Format mesh reader.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Hex_Format.hh"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

namespace rtt_meshReaders
{

Hex_Format::Hex_Format(std::string filename)
    : meshfile_name(filename)
{

    // Open the mesh file for read.
    std::ifstream meshfile(filename.c_str(), std::ios::in);
    if (!meshfile)
	Insist(false,"Could not open mesh-file!");

    // Check to make sure this is a CIC-19 Hex format file.
    std::string chdum;
    meshfile >> chdum;
    if ( chdum != keyword() )
	Insist(false,"Not a CIC-19 Hex Mesh File!");

    // Read in the dimensions of the problem.
    meshfile >> npoints >> ncells >> nvrtx >> nvrpf >> ndim >> nvb_faces
	     >> nrb_faces >> nmat;

    std::cout << "npoints= " << npoints << std::endl;
    std::cout << "ncells= " << ncells << std::endl;
    std::cout << "nvrtx= " << nvrtx << std::endl;
    std::cout << "nvrpf= " << nvrpf << std::endl;
    std::cout << "ndim= " << ndim  << std::endl;
    std::cout << "nvb_faces= " << nvb_faces << std::endl;
    std::cout << "nrb_faces= " <<  nrb_faces << std::endl;
    std::cout << "nmat= " <<  nmat << std::endl;

    // Read the point coordinates data.
    point_coords.resize(npoints);
    std::cout << "Capacity of Point Coords: " << point_coords.capacity() << std::endl;
    std::cout << "Size of Point Coords: " << point_coords.size() << std::endl;
    for (int i=0; i<npoints; i++)
    {
	point_coords[i].resize(ndim);
	if (ndim == 1)
	    meshfile >> point_coords[i][0];
	else if (ndim == 2)
	    meshfile >> point_coords[i][0] >> point_coords[i][1];
	else if (ndim == 3)
	    meshfile >> point_coords[i][0] >> point_coords[i][1]
		     >> point_coords[i][2];
	else
	    Insist(false,"Dimension index out of range!");  
    }

    std::cout << "Read coordinates." << std::endl;

    // Read in the mesh connectivity.
    ipar.resize(ncells);
    for (int i=0; i<ncells; i++)
    {
	ipar[i].resize(nvrtx);
	if (ndim == 1) 
	    meshfile >> ipar[i][0] >> ipar[i][1];
	else if (ndim ==2)
	    meshfile >> ipar[i][0] >> ipar[i][1] >> ipar[i][2] >> ipar[i][3];
	else if (ndim ==3)
	    meshfile >> ipar[i][0] >> ipar[i][1] >> ipar[i][2] >> ipar[i][3]
		     >> ipar[i][4] >> ipar[i][5] >> ipar[i][6] >> ipar[i][7];
	else
	    Insist(false,"Dimension index out of range!");
	for (int j=0; j<nvrtx; j++)
	    ipar[i][j] = ipar[i][j]-1;
    }
    std::cout << "Read elements." << std::endl;

    // Read in the mesh interior-region data.
    imat_index.resize(ncells);
    for (int i=0; i<ncells; i++)
	meshfile >> imat_index[i];
    std::cout << "Read material flag." << std::endl;

    // Read in the mesh vacuum boundary data.
    ipar_vb.resize(nvb_faces);
    irgn_vb_index.resize(nvb_faces);
    for (int i=0; i<nvb_faces; i++)
    {
	ipar_vb[i].resize(nvrpf);
	if (ndim == 1)
	    meshfile >> ipar_vb[i][0] >> irgn_vb_index[i];
	else if (ndim ==2)
	    meshfile >> ipar_vb[i][0] >> ipar_vb[i][1] >> irgn_vb_index[i];
	else if (ndim ==3)
	    meshfile >> ipar_vb[i][0] >> ipar_vb[i][1] >>ipar_vb[i][2] >>
		ipar_vb[i][3] >> irgn_vb_index[i];
	else
	    Insist(false,"Dimension index out of range!");
	for (int j=0; j<nvrpf; j++) 
	    ipar_vb[i][j] = ipar_vb[i][j]-1;
    }
    std::cout << "Read vacuum boundary data." << std::endl;

    // Read in the mesh reflective boundary data.
    ipar_rb.resize(nrb_faces);
    for (int i=0; i<nrb_faces; i++)
    {
	ipar_rb[i].resize(nvrpf);
	if (ndim ==1)
	    meshfile >> ipar_rb[i][0];
	else if(ndim ==2)
	    meshfile >> ipar_rb[i][0] >> ipar_rb[i][1];
	else if (ndim == 3)
	    meshfile >> ipar_rb[i][0] >> ipar_rb[i][1] >>ipar_rb[i][2] >>
		ipar_vb[i][3];
	else
	    Insist(false,"Dimension index out of range!");
	for (int j=0; j<nvrpf; j++) 
	    ipar_rb[i][j] = ipar_rb[i][j]-1;
    }
    std::cout << "Read reflective boundary data." << std::endl;

    Ensure (invariant());
}

std::vector<std::vector<int> > Hex_Format::get_element_nodes() const
{
    // Collate the interior, vacuum, and reflective mesh elements into
    // one vector. Note that the order is important as we will rely
    // on it later to output the element set data.

    // Alternatively, the private data of the class could be
    // changed so that the work done here is done in the constructor.
    // This would be more efficient if this is going to be
    // used repetively.
    std::vector<std::vector<int> > result;
    result.resize(ncells+nvb_faces+nrb_faces);
    for (int i=0; i<ncells; i++)
	result.push_back(ipar[i]);
    for (int i=0; i<nvb_faces; i++)
	result.push_back(ipar_vb[i]);
    for (int i=0; i<nrb_faces; i++)
	result.push_back(ipar_rb[i]);
    return result;
}

std::vector<Element_Definition::Element_Type> Hex_Format::get_element_types() const
{
    Element_Definition::Element_Type d1, d2;
    switch ( ndim )
    { 
    case (1) :
	d1 = Element_Definition::BAR_2;
	d2 = Element_Definition::NODE;
	break;
    case (2) :
	d1 = Element_Definition::QUAD_4;
	d2 = Element_Definition::BAR_2;
	break;
    case (3) :
	d1 = Element_Definition::HEXA_8;
	d2 = Element_Definition::QUAD_4;
	break;
    default :
	Insist(false,"Dimension index out of range!");
    }
    std::vector<Element_Definition::Element_Type> tmp;
    for(int i=0; i<ncells; i++)
	tmp.push_back(d1);
    for(int i=0; i<nvb_faces+nrb_faces; i++)
	tmp.push_back(d2);
    return tmp;
}

std::map<std::string, std::set<int> > Hex_Format::get_element_sets() const
{
    // Alternatively, the private data of the class could be
    // changed so that the work done here is done in the constructor.
    // This would be more efficient if this is going to be
    // used repetively.
    typedef std::map<std::string, std::set<int> > resultT;
    resultT result;
    std::vector<int> tmp;
    std::set<int> rgn_index;
    std::ostringstream os_chdum;
    std::set<int> stmp;

    // Create a vacuum boundary set. Note that this depends on
    // the elements being stored in a specific order.
    for (int i=ncells; i<ncells+nvb_faces; i++)
	stmp.insert(i);
    result.insert(resultT::value_type("Vacuum_Boundary", stmp));
    
    // Create a reflective boundary set. Note that this depends on
    // the elements being stored in a specific order.
    stmp.clear();
    for (int i=ncells+nvb_faces; i<ncells+nvb_faces+nrb_faces; i++)
	stmp.insert(i);
    result.insert(resultT::value_type("Reflective_Boundary", stmp));

    // Create sets for all the interior mesh regions.
    // This loops over the whole mesh number_of_mesh_regions times. Could
    // be made to do it more efficiently in one loop? Note that this 
    // depends on  the elements being stored in a specific order.
    rgn_index = std::set<int>(imat_index.begin(), imat_index.end());
    for (std::set<int>::iterator i=rgn_index.begin(); i != rgn_index.end(); i++)
    {
	os_chdum.clear();
	os_chdum << "Interior_Region_" << *i;
	stmp.clear();
	for (int j=0; j<ncells; j++)
	    if (imat_index[j] == *i) stmp.insert(j);
	result.insert(resultT::value_type(os_chdum.str(), stmp));
    }	

    // Create sets for all the vacuum boundary regions.
    // This loops over the whole mesh number_of_vb_regions times. Could
    // be made to do it more efficiently in one loop? Note that 
    // this depends on the elements being stored in a specific order.
    rgn_index = std::set<int>(irgn_vb_index.begin(), irgn_vb_index.end());
    for (std::set<int>::iterator i=rgn_index.begin(); i != rgn_index.end(); i++)
    {
	os_chdum.clear();	
	os_chdum << "Vacuum_Boundary_Region_" << *i;
	stmp.clear();
	for (int j=0; j<ncells; j++)
	    if (irgn_vb_index[j] == *i) stmp.insert(j);
	result.insert(resultT::value_type( os_chdum.str(), stmp));
    }	
    return result;
}

bool Hex_Format::invariant() const
{
    bool ldum =
	(npoints > 0) && (ncells > 0) 
	&& ( (ndim == 3 && nvrtx == 8 && nvrpf == 4) || 
	     (ndim == 2 && nvrtx == 4 && nvrpf == 2) ||
	     (ndim == 1 && nvrtx == 2 && nvrpf == 1) )
	&& (nvb_faces >= 0) && (nrb_faces >= 0)
	&& (nmat > 0) 
	&& (point_coords.size() == npoints)
	&& (ipar.size() == ncells)
	&& (imat_index.size() == ncells)
	&& (irgn_vb_index.size() == nvb_faces)
	&& (ipar_vb.size() == nvb_faces)
	&& (ipar_rb.size() == nrb_faces);
    return ldum;
}

} // end namespace rtt_meshReaders


//---------------------------------------------------------------------------//
//                              end of Hex_Format.cc
//---------------------------------------------------------------------------//
