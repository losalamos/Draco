//----------------------------------*-C++-*--------------------------------//
// RTT_Mesh_Reader.cc
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/RTT_Mesh_Reader.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for RTT_Mesh_Reader library.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "RTT_Mesh_Reader.hh"

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Transforms the RTT_Format data to the CYGNUS format.
 */
void RTT_Mesh_Reader::transform2CYGNUS()
{
    rtt_meshReaders::Element_Definition::Element_Type cell_defs;
    std::vector<rtt_meshReaders::Element_Definition::Element_Type> cell_types;
    rtt_dsxx::SP<rtt_meshReaders::Element_Definition> cell;
    vector_int new_side_types;
    vector_vector_int new_ordered_sides;
    vector_vector_int cell_side_types(rttMesh->get_dims_ncell_defs());
    vector_vector_vector_int 
        cell_ordered_sides(rttMesh->get_dims_ncell_defs());

    for (int cd = 0; cd < rttMesh->get_dims_ncell_defs(); cd++)
    {
        string cell_name = rttMesh->get_cell_defs_name(cd);

	if (cell_name == "point")
	    cell_defs = rtt_meshReaders::Element_Definition::NODE;
	else if (cell_name == "line")
	    cell_defs = rtt_meshReaders::Element_Definition::BAR_2;
	else if (cell_name == "triangle")
	    cell_defs = rtt_meshReaders::Element_Definition::TRI_3;
	else if (cell_name == "quad")
	    cell_defs = rtt_meshReaders::Element_Definition::QUAD_4;
	else if (cell_name == "tetrahedron")
	    cell_defs = rtt_meshReaders::Element_Definition::TETRA_4;
	else if (cell_name == "quad_pyr")
	    cell_defs = rtt_meshReaders::Element_Definition::PYRA_5;
	else if (cell_name == "tri_prism")
	    cell_defs = rtt_meshReaders::Element_Definition::PENTA_6;
	else if (cell_name == "hexahedron")
	    cell_defs = rtt_meshReaders::Element_Definition::HEXA_8;
	else
	    throw std::runtime_error("Unrecognized cell definition");

	cell_types.push_back(cell_defs);
	cell = new rtt_meshReaders::Element_Definition(cell_defs);
	new_side_types.resize(cell->get_number_of_sides());
	new_ordered_sides.resize(cell->get_number_of_sides());
	for (int s = 0; s < cell->get_number_of_sides(); s++)
	{
	    new_side_types[s] = (std::find(cell_types.begin(), 
					   cell_types.end(),
					   cell->get_side_type(s).get_type())
				 - cell_types.begin());
	    new_ordered_sides[s] = cell->get_side_nodes(s);
	}
	cell_side_types[cd] = new_side_types;
	cell_ordered_sides[cd] = new_ordered_sides;
    }
    rttMesh->reformatData(cell_side_types, cell_ordered_sides);

    // Load the element types vector.
    for (int s = 0; s < rttMesh->get_dims_nsides(); s++)
        element_types.push_back(cell_types[rttMesh->get_sides_type(s)]);
    for (int c = 0; c < rttMesh->get_dims_ncells(); c++)
        element_types.push_back(cell_types[rttMesh->get_cells_type(c)]);
}
/*!
 * \brief Returns the node numbers associated with each element (i.e., sides
 *        and cells).
 * \return The node numbers.
 */
std::vector<std::vector<int> > RTT_Mesh_Reader::get_element_nodes() const
{
    vector_vector_int element_nodes(rttMesh->get_dims_nsides() + 
				    rttMesh->get_dims_ncells());

    for (int i = 0; i < rttMesh->get_dims_nsides(); i++)
	element_nodes[i] = rttMesh->get_sides_nodes(i);

    int nsides = rttMesh->get_dims_nsides();
    for (int i = 0; i < rttMesh->get_dims_ncells(); i++)
	element_nodes[i + nsides] = rttMesh->get_cells_nodes(i);

    return element_nodes;
}
/*!
 * \brief Returns the nodes associated with each node_flag_type_name and
 *        node_flag_name combination.
 * \return The nodes associated with each node_flag_type_name/node_flag_name 
 *         combination.
 */
std::map<std::string, std::set<int> > RTT_Mesh_Reader::get_node_sets() const
{
    std::map<string, set_int > node_sets;
    string flag_types_and_names;

    // loop over the number of node flag types.
    for (int type = 0; type < rttMesh->get_dims_nnode_flag_types(); type++)
    {
        // loop over the number of node flags for this type.
        for (int flag = 0; flag < rttMesh->get_dims_nnode_flags(type); flag++)
	{
            set_int node_flags;
            flag_types_and_names =  rttMesh->get_node_flags_flag_type(type);
	    flag_types_and_names.append("/");
	    flag_types_and_names += rttMesh->get_node_flags_flag_name(type, 
								      flag);
	    int flag_number = rttMesh->get_node_flags_flag_number(type, flag);
            // loop over the nodes.
	    for (int node = 0; node < rttMesh->get_dims_nnodes(); node++)
	    {
	        if (flag_number == rttMesh->get_nodes_flags(node,type))
		   node_flags.insert(node);
	    }
	    node_sets.insert(std::make_pair(flag_types_and_names, node_flags));
	}
    }
    return node_sets;
}
/*!
 * \brief Returns the elements (i.e., sides and cells) associated with each 
 *        flag_type_name and flag_name combination for the sides and cells
 *        read from the mesh file data.
 * \return The elements associated with each flag_type_name/flag_name 
 *         combination.
 */
std::map<std::string, std::set<int> > RTT_Mesh_Reader::get_element_sets() const
{
    std::map<string, set_int > element_sets;
    string flag_types_and_names;

    // loop over the number of side flag types.
    for (int type = 0; type < rttMesh->get_dims_nside_flag_types(); type++)
    {
        // loop over the number of side flags for this type.
        for (int flag = 0; flag < rttMesh->get_dims_nside_flags(type); flag++)
	{
            set_int side_flags;
            flag_types_and_names = rttMesh->get_side_flags_flag_type(type);
	    flag_types_and_names.append("/");
	    flag_types_and_names +=rttMesh->get_side_flags_flag_name(type, 
								     flag);
	    int flag_number =rttMesh->get_side_flags_flag_number(type, flag);
            // loop over the sides.
	    for (int side = 0; side < rttMesh->get_dims_nsides(); side++)
	    {
	        if (flag_number == rttMesh->get_sides_flags(side, type))
		   side_flags.insert(side);
	    }
	    element_sets.insert(std::make_pair(flag_types_and_names,
					       side_flags));
	}
        
    }
    int nsides = rttMesh->get_dims_nsides();
    // loop over the number of cell flag types.
    for (int type = 0; type < rttMesh->get_dims_ncell_flag_types(); type++)
    {
        // loop over the number of cell flags for this type.
        for (int flag = 0; flag < rttMesh->get_dims_ncell_flags(type); flag++)
	{
            set_int cell_flags;
            flag_types_and_names =  rttMesh->get_cell_flags_flag_type(type);
	    flag_types_and_names.append("/");
	    flag_types_and_names += rttMesh->get_cell_flags_flag_name(type, 
								      flag);
	    int flag_number = rttMesh->get_cell_flags_flag_number(type, flag);
            // loop over the cells.
	    for (int cell = 0; cell < rttMesh->get_dims_ncells(); cell++)
	    {
	        if (flag_number == rttMesh->get_cells_flags(cell, type))
		   cell_flags.insert(cell + nsides);
	    }
	    // Allow the possibility that the cells could haved identical 
	    // flags as the sides.
	    if (element_sets.count(flag_types_and_names) != 0)
	    {
	        set_int side_set =
		    element_sets.find(flag_types_and_names)->second;
		for (set_int::const_iterator side_set_itr = side_set.begin();
		     side_set_itr != side_set.end(); side_set_itr++)
		    cell_flags.insert(* side_set_itr);
		element_sets.erase(flag_types_and_names);
	    }
	    element_sets.insert(std::make_pair(flag_types_and_names,
					       cell_flags));
	}
    }

    return element_sets;
}
/*!
 * \brief Performs a basic sanity check on the mesh file data.
 * \return Acceptablity of the mesh file data.
 */
bool RTT_Mesh_Reader::invariant() const
{
    bool test =  (rttMesh->get_dims_ndim() > 0) && 
                 (rttMesh->get_dims_nnodes() > 0) && 
                 (rttMesh->get_dims_nsides() > 0) && 
                 (rttMesh->get_dims_ncells() > 0) &&
                 (rttMesh->get_dims_ncell_defs() > 0) && 
                 (rttMesh->get_dims_nside_types() > 0) &&
                 (rttMesh->get_dims_ncell_types() > 0) &&
                 (rttMesh->get_dims_ncell_defs() >= 
		     rttMesh->get_dims_nside_types()) &&
		 (rttMesh->get_dims_ncell_defs() >= 
		     rttMesh->get_dims_ncell_types());
    return test;
}

} // end namespace rtt_RTT_Format_Reader

//---------------------------------------------------------------------------//
//                   end of RTT_Format_Reader/RTT_Mesh_Reader.cc
//---------------------------------------------------------------------------//
