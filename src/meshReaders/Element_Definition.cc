//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/Element_Definition.cc
 * \author John McGhee
 * \date   Fri Feb 25 10:03:18 2000
 * \brief  Provides some descriptive information for the
 *         standard mesh elements.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Element_Definition.hh"
#include "ds++/Assert.hh"

namespace rtt_meshReaders
{

Element_Definition::Element_Definition(const Element_Type &type_)

    : type(type_)
{
    switch( type )
    {

    case NODE :
	construct_node();
	break;
	
    case BAR_2 :
    case BAR_3 :
	construct_bar();
	break;

    case TRI_3 :
    case TRI_6 :
	construct_tri();
	break;

    case QUAD_4 :
    case QUAD_8 :
    case QUAD_9 :
	construct_quad();
	break;

    case TETRA_4  :
    case TETRA_10 :
	construct_tetra();
	break;

    case PYRA_5  :
    case PYRA_14 :
	construct_pyra();
	break;

    case PENTA_6  :
    case PENTA_15 :
    case PENTA_18 :
	construct_penta();
	break;	

    case HEXA_8  :
    case HEXA_20 :
    case HEXA_27 :
	construct_hexa();
	break;
	
    default :
	
	throw std::runtime_error("#1 Unrecognized Element-Type Flag");
	
    }
    
    Ensure( invariant_satisfied() );
}

bool Element_Definition::invariant_satisfied() const
{
    
    bool ldum = name.empty() == false;

    if (type == NODE)
    {
	ldum = ldum && 
	    dimension         == 0 && number_of_nodes  == 1 &&
	    elem_defs.size()  == 0 && side_type.size() == 0 &&
	    side_nodes.size() == 0 && number_of_sides  == 0 ;
    }
    else
    {
	ldum = ldum &&
	    dimension > 0 && 
	    number_of_nodes > dimension &&
	    elem_defs.size() > 0 &&
	    node_loc.size() == number_of_nodes;
	for (int i=0; i < elem_defs.size(); i++)
	{
	    ldum = ldum && elem_defs[i].dimension == dimension-1;
	}
	ldum = ldum && 
	    number_of_sides <= number_of_nodes &&
	    number_of_sides > dimension; 
	ldum == ldum && side_type.size() == number_of_sides;
	ldum == ldum && side_nodes.size() == number_of_sides;
        for (int i=0; i < side_nodes.size(); i++)
	{
	    ldum = ldum && side_nodes[i].size() > 0;
	    ldum = ldum && side_nodes[i].size() == 
		elem_defs[ side_type[i] ].number_of_nodes;
	    for (int j=0; j < side_nodes[i].size(); j++)
	    {
		ldum = ldum && side_nodes[i][j] >= 0 &&
		    side_nodes[i][j] < number_of_nodes;
		ldum = ldum && node_loc[ side_nodes[i][j] ] ==
		    elem_defs[ side_type[i] ].node_loc[j];
	    }
	}
    }
    return ldum;
}

void Element_Definition::construct_node()
{
    name = "NODE";
    dimension=0;
    number_of_sides=0;
    number_of_nodes=1;
    node_loc.push_back(CORNER);
}

void Element_Definition::construct_bar()
{
    std::vector<int> tmp;
    dimension=1;
    number_of_sides=2;
    tmp.clear();
    tmp.push_back(0);
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    side_nodes.push_back(tmp);
    for (int i=0; i < 2; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    {
    case BAR_2 :
	name = "BAR_2";
	number_of_nodes=2;
	break;
    case BAR_3 :
	name = "BAR_3";
	number_of_nodes=3;
	node_loc.push_back(EDGE);
	break;
    default :
	throw std::runtime_error("#2 Unrecognized Element-Type Flag");
    }
    elem_defs.push_back(Element_Definition(NODE));
    for (int i = 0; i < number_of_sides; i++)
	side_type.push_back(0);
}

void Element_Definition::construct_tri()
{
    std::vector<int> tmp;
    dimension=2;
    number_of_sides=3;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(1); 
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 0;
    side_nodes.push_back(tmp);
    for (int i=0; i < 3; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    {
    case  TRI_3 :
	name = "TRI_3";
	elem_defs.push_back(Element_Definition(BAR_2));
	number_of_nodes=3;
	break;
    case TRI_6 :
	name = "TRI_6";
	elem_defs.push_back(Element_Definition(BAR_3));
	number_of_nodes=6;
	side_nodes[0].push_back(3);
	side_nodes[1].push_back(4);
	side_nodes[2].push_back(5);
	for (int i=0; i < 3; i++)
	    node_loc.push_back(EDGE);
	break;
    default :
	throw std::runtime_error("#3 Unrecognized Element-Type Flag");
    } 
    for (int i = 0; i < number_of_sides; i++)
	side_type.push_back(0);
}

void Element_Definition::construct_quad()
{
    std::vector<int> tmp;
    dimension=2;
    number_of_sides=4;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(1);    
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 3;
    side_nodes.push_back(tmp);
    tmp[0] = 3;
    tmp[1] = 0;
    side_nodes.push_back(tmp);
    for (int i=0; i < 4; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    {
    case QUAD_4 :
	name = "QUAD_4";
	number_of_nodes=4;
	elem_defs.push_back(Element_Definition(BAR_2));
	break;
    case QUAD_8 :
    case QUAD_9 :
	elem_defs.push_back(Element_Definition(BAR_3));
	for (int i=0; i < 4; i++)
	    node_loc.push_back(EDGE);
	side_nodes[0].push_back(4);
	side_nodes[1].push_back(5);
	side_nodes[2].push_back(6);
	side_nodes[3].push_back(7);
	switch ( type )
	{
	case QUAD_8 :
	    name = "QUAD_8";
	    number_of_nodes=8;
	    break;
	case QUAD_9 :
	    name = "QUAD_9";
	    number_of_nodes=9;
	    node_loc.push_back(FACE);
	    break;
	default :
	    throw std::runtime_error("#4 Unrecognized Element-Type Flag");
	}
	break;
    default :
	throw std::runtime_error("#5 Unrecognized Element-Type Flag");
    }
    for (int i = 0; i < number_of_sides; i++)
	side_type.push_back(0);
}

void Element_Definition::construct_tetra()
{
    std::vector<int> tmp;
    dimension=3;
    number_of_sides=4;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(2); 
    tmp.push_back(1);    
    side_nodes.push_back(tmp);
    tmp[0] = 0;
    tmp[1] = 1;
    tmp[2] = 3;
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    tmp[2] = 3;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 0;
    tmp[2] = 3;
    side_nodes.push_back(tmp);
    for (int i=0; i < 4; i++)
	node_loc.push_back(CORNER);	
    switch ( type )
    {
    case TETRA_4 :
	name = "TETRA_4";
	number_of_nodes=4;
	elem_defs.push_back(Element_Definition(TRI_3));
	break;
    case TETRA_10 :
	name = "TETRA_10";
	number_of_nodes=10;
	elem_defs.push_back(Element_Definition(TRI_6));
	
	side_nodes[0].push_back(6);
	side_nodes[0].push_back(5);
	side_nodes[0].push_back(4);
	
	side_nodes[1].push_back(4);
	side_nodes[1].push_back(8);
	side_nodes[1].push_back(7);
	
	side_nodes[2].push_back(5);
	side_nodes[2].push_back(9);
	side_nodes[2].push_back(8);
	
	side_nodes[3].push_back(6);
	side_nodes[3].push_back(7);
	side_nodes[3].push_back(9);
	
	for (int i=0; i < 6; i++)
	    node_loc.push_back(EDGE);
	break;
    default :
	throw std::runtime_error("#6 Unrecognized Element-Type Flag");
    }
    for (int i = 0; i < number_of_sides; i++)
	side_type.push_back(0);
}

void Element_Definition::construct_pyra()
{
    std::vector<int> tmp;
    dimension=3;
    number_of_sides=5;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(3); 
    tmp.push_back(2);
    tmp.push_back(1);
    side_nodes.push_back(tmp);
    tmp.pop_back();
    tmp[0] = 0;
    tmp[1] = 1;
    tmp[2] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    tmp[2] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 3;
    tmp[2] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 3;
    tmp[1] = 0;
    tmp[2] = 4;
    side_nodes.push_back(tmp);
    for (int i=0; i < 5; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    {
    case PYRA_5 :
	name = "PYRA_5";
	number_of_nodes=5;
	elem_defs.push_back(Element_Definition(QUAD_4));
	elem_defs.push_back(Element_Definition(TRI_3));
	break;
    case PYRA_14 :
	name = "PYRA_14";
	number_of_nodes=14;
	elem_defs.push_back(Element_Definition(QUAD_8));
	elem_defs.push_back(Element_Definition(TRI_6));
	for (int i=0; i < 8; i++)
	    node_loc.push_back(EDGE);
	node_loc.push_back(CELL);
	
	side_nodes[0].push_back(8);
	side_nodes[0].push_back(7);
	side_nodes[0].push_back(6);
	side_nodes[0].push_back(5);
	
	side_nodes[1].push_back(5);
	side_nodes[1].push_back(10);
	side_nodes[1].push_back(9);
	
	side_nodes[2].push_back(6);
	side_nodes[2].push_back(11);
	side_nodes[2].push_back(10);
	
	side_nodes[3].push_back(7);
	side_nodes[3].push_back(12);
	side_nodes[3].push_back(11);
	
	side_nodes[4].push_back(8);
	side_nodes[4].push_back(9);
	side_nodes[4].push_back(12);
	break;
    default :
	throw std::runtime_error("#7 Unrecognized Element-Type Flag");
    }
    side_type.push_back(0);
    for (int i = 1; i < number_of_sides; i++)
	side_type.push_back(1);
}

void Element_Definition::construct_penta()
{
    std::vector<int> tmp;
    dimension=3;
    number_of_sides=5;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(1); 
    tmp.push_back(4); 
    tmp.push_back(3);   
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    tmp[2] = 5;
    tmp[3] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 0;
    tmp[2] = 3;
    tmp[3] = 5;
    side_nodes.push_back(tmp);
    tmp.pop_back();
    tmp[0] = 0;
    tmp[1] = 2;
    tmp[2] = 1;
    side_nodes.push_back(tmp);
    tmp[0] = 3;
    tmp[1] = 4;
    tmp[2] = 5;
    side_nodes.push_back(tmp);
    for (int i=0; i < 6; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    { 
    case PENTA_6 :
	name = "PENTA_6";
	number_of_nodes=6;
	elem_defs.push_back(Element_Definition(QUAD_4));
	elem_defs.push_back(Element_Definition(TRI_3));
	break;
    case PENTA_15 :
    case PENTA_18 :
	for (int i=0; i < 9; i++)
	    node_loc.push_back(EDGE);
	
	side_nodes[0].push_back(6);
	side_nodes[0].push_back(10);
	side_nodes[0].push_back(12);
	side_nodes[0].push_back(9);
	
	side_nodes[1].push_back(7);
	side_nodes[1].push_back(11);
	side_nodes[1].push_back(13);
	side_nodes[1].push_back(10);
	
	side_nodes[2].push_back(8);
	side_nodes[2].push_back(9);
	side_nodes[2].push_back(14);
	side_nodes[2].push_back(11);
	
	side_nodes[3].push_back(8);
	side_nodes[3].push_back(7);
	side_nodes[3].push_back(6);
	
	side_nodes[4].push_back(12);
	side_nodes[4].push_back(13);
	side_nodes[4].push_back(14);
	switch ( type )
	{
	case PENTA_15 :
	    name = "PENTA_15";
	    number_of_nodes=15;
	    elem_defs.push_back(Element_Definition(QUAD_8));
	    elem_defs.push_back(Element_Definition(TRI_6));
	    break;
	case PENTA_18 :
	    name = "PENTA_18";
	    for (int i=0; i < 3; i++)
		node_loc.push_back(FACE);
	    number_of_nodes=18;
	    elem_defs.push_back(Element_Definition(QUAD_9));
	    elem_defs.push_back(Element_Definition(TRI_6));
	    side_nodes[0].push_back(15);
	    side_nodes[1].push_back(16);
	    side_nodes[2].push_back(17);
	    break;
	default :
	    throw std::runtime_error("#8 Unrecognized Element-Type Flag");
	}
	break;
    default :
	throw std::runtime_error("#9 Unrecognized Element-Type Flag");
    }
    for (int i = 0; i < 3; i++)
	side_type.push_back(0);
    for (int i = 3; i < number_of_sides; i++)
	side_type.push_back(1);
}

void Element_Definition::construct_hexa()
{
    std::vector<int> tmp;
    dimension=3;
    number_of_sides=6;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(3); 
    tmp.push_back(2); 
    tmp.push_back(1);   
    side_nodes.push_back(tmp);
    tmp[0] = 0;
    tmp[1] = 1;
    tmp[2] = 5;
    tmp[3] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    tmp[2] = 6;
    tmp[3] = 5;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 3;
    tmp[2] = 7;
    tmp[3] = 6;
    side_nodes.push_back(tmp);
    tmp[0] = 0;
    tmp[1] = 4;
    tmp[2] = 7;
    tmp[3] = 3;
    side_nodes.push_back(tmp);
    tmp[0] = 4;
    tmp[1] = 5;
    tmp[2] = 6;
    tmp[3] = 7;
    side_nodes.push_back(tmp);
    for (int i=0; i < 8; i++)
	node_loc.push_back(CORNER);
    switch ( type )
    {
    case HEXA_8 :
	name = "HEXA_8";
	number_of_nodes=8;
	elem_defs.push_back(Element_Definition(QUAD_4));
	break;
    case HEXA_20 :
    case HEXA_27 :
	for (int i=0; i < 12; i++)
	    node_loc.push_back(EDGE);
	
	side_nodes[0].push_back(11);
	side_nodes[0].push_back(10);
	side_nodes[0].push_back(9);
	side_nodes[0].push_back(8);
	
	side_nodes[1].push_back(8);
	side_nodes[1].push_back(13);
	side_nodes[1].push_back(16);
	side_nodes[1].push_back(12);
	
	side_nodes[2].push_back(9);
	side_nodes[2].push_back(14);
	side_nodes[2].push_back(17);
	side_nodes[2].push_back(13);
	
	side_nodes[3].push_back(10);
	side_nodes[3].push_back(15);
	side_nodes[3].push_back(18);
	side_nodes[3].push_back(14);
	
	side_nodes[4].push_back(12);
	side_nodes[4].push_back(19);
	side_nodes[4].push_back(15);
	side_nodes[4].push_back(11);
	
	side_nodes[5].push_back(16);
	side_nodes[5].push_back(17);
	side_nodes[5].push_back(18);
	side_nodes[5].push_back(19);
	switch ( type )
	{
	case HEXA_20 :
	    name = "HEXA_20";
	    number_of_nodes=20;
	    elem_defs.push_back(Element_Definition(QUAD_8));
	    break;
	case HEXA_27 :
	    name = "HEXA_27";
	    for (int i=0; i < 6; i++)
		node_loc.push_back(FACE);
	    node_loc.push_back(CELL);
	    number_of_nodes=27;
	    elem_defs.push_back(Element_Definition(QUAD_9));
	    for (int i=0; i < number_of_sides; i++)
		side_nodes[i].push_back(i+20);
	    break;
	default :
	    throw std::runtime_error("#10 Unrecognized Element-Type Flag");
	}
	break;
    default :
	throw std::runtime_error("#11 Unrecognized Element-Type Flag");
    }
    for (int i = 0; i < number_of_sides; i++)
	side_type.push_back(0);	
}

} // end namespace rtt_meshReaders

//---------------------------------------------------------------------------//
//                              end of Element_Definition.cc
//---------------------------------------------------------------------------//
