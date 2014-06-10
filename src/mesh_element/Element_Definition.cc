//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_element/Element_Definition.cc
 * \author John McGhee
 * \date   Fri Feb 25 10:03:18 2000
 * \brief  Provides some descriptive information for the
 *         standard mesh elements.
 * \note   Copyright (C) 2000-2014 Los Alamos National Security, LLC. 
 *         All rights reserved. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "Element_Definition.hh"

namespace rtt_mesh_element
{

//---------------------------------------------------------------------------//
Element_Definition::Element_Definition( Element_Type const & type_ )
    : name(),
      type(type_),
      dimension(0),
      number_of_nodes(0),
      number_of_sides(0),
      elem_defs(),
      side_type(),
      side_nodes(),
      node_loc()
{
    switch( type )
    {
	
    case NODE :
	construct_node();
	break;
	
    case BAR_2 :
    case BAR_3 :
	construct_bar();
        Ensure( invariant_satisfied() );
	break;

    case TRI_3 :
    case TRI_6 :
	construct_tri();
        Ensure( invariant_satisfied() );
	break;

    case QUAD_4 :
    case QUAD_8 :
    case QUAD_9 :
	construct_quad();
        Ensure( invariant_satisfied() );
	break;

    case TETRA_4  :
    case TETRA_10 :
	construct_tetra();
        Ensure( invariant_satisfied() );
	break;

    case PYRA_5  :
    case PYRA_14 :
	construct_pyra();
        Ensure( invariant_satisfied() );
	break;

    case PENTA_6  :
    case PENTA_15 :
    case PENTA_18 :
	construct_penta();
        Ensure( invariant_satisfied() );
	break;	

    case HEXA_8  :
    case HEXA_20 :
    case HEXA_27 :
	construct_hexa();
        Ensure( invariant_satisfied() );
	break;	

    case POLYGON :
        dimension = 2;
	break;	

    default :
	Insist(false,"Unrecognized Element-Type Flag");
	
    }

//     number_of_face_nodes.resize(side_type.size());
//     for (unsigned s=0; s < side_type.size(); ++s)
//         number_of_face_nodes[s] = get_side_type(s).get_number_of_nodes();

//     face_nodes.resize(side_nodes.size());

//     for (unsigned s=0; s < face_nodes.size(); ++s)
//     {
//         std::vector<int> nodes(get_side_nodes(s)); 
//         face_nodes[s].resize(nodes.size());

//         for (unsigned n=0; n < nodes.size(); ++n)
//             face_nodes[s][n] = nodes[n];
//     }
    
}
    
//---------------------------------------------------------------------------//
Element_Definition::
Element_Definition( std::string  name_,
                    size_t dimension_,
                    size_t number_of_nodes_,
                    size_t number_of_sides_,
                    std::vector< Element_Definition >  const & elem_defs_,
                    std::vector< int >                 const & side_type_,
                    std::vector< std::vector<size_t> > const & side_nodes_,
                    std::vector< Node_Location >       const & node_loc_ )
    :
    name(name_),
    type(POLYGON),
    dimension(dimension_),
    number_of_nodes(number_of_nodes_),
    number_of_sides(number_of_sides_),
    elem_defs(elem_defs_),
    side_type(side_type_),
    side_nodes(side_nodes_),
    node_loc(node_loc_)
{
    Require( number_of_nodes_>0 );
    
#if DBC & 1
    for (unsigned i=0; i<elem_defs_.size(); ++i)
    {
        Require( elem_defs_[i].get_dimension()+1==dimension_ );
    }
#endif
    Require( side_type_.size()==number_of_sides_ );
#if DBC & 1
    for (unsigned i=0; i<number_of_sides_; ++i)
    {
        Require( static_cast<unsigned>(side_type_[i])<elem_defs.size() );
    }
#endif
    Require( side_nodes_.size()==number_of_sides_ );
    Require( node_loc_.size()==number_of_nodes_ );
#if DBC & 1
    for (unsigned i=0; i<number_of_sides_; ++i)
    {
        Require( side_nodes_[i].size() ==
                 elem_defs_[side_type_[i]].get_number_of_nodes() );
        
        for (unsigned j=0; j<side_nodes_[i].size(); ++j)
        {
            Require( static_cast<unsigned>(side_nodes_[i][j])
                     < number_of_nodes_ );
            
            Require( elem_defs_[side_type_[i]].get_node_location(j) ==
                     node_loc_[side_nodes_[i][j]] );
        }
    }
#endif

     Ensure( invariant_satisfied() );

     Ensure( get_type()==Element_Definition::POLYGON );
     Ensure( get_name()==name_  );
     Ensure( get_dimension()==dimension_  );
     Ensure( get_number_of_nodes()==number_of_nodes_  );
     Ensure( get_number_of_sides()==number_of_sides_  );
#if DBC & 4
     for (unsigned i=0; i<number_of_nodes; ++i)
     {
         Ensure( get_node_location(i)==node_loc_[i]  );
     }
     for (unsigned i=0; i<number_of_sides; ++i)
     {
//         Ensure( get_side_type(i)==elem_defs_[side_type_[i]]  );
         Ensure( get_side_nodes(i)==side_nodes_[i]  );
     }
#endif

//     number_of_face_nodes.resize(side_type.size());
//     for (unsigned s=0; s < side_type.size(); ++s)
//         number_of_face_nodes[s] = get_side_type(s).get_number_of_nodes();

//     face_nodes.resize(side_nodes.size());

//     for (unsigned s=0; s < face_nodes.size(); ++s)
//     {
//         std::vector<int> nodes(get_side_nodes(s)); 
//         face_nodes[s].resize(nodes.size());

//         for (unsigned n=0; n < nodes.size(); ++n)
//             face_nodes[s][n] = nodes[n];
//     }
}

//---------------------------------------------------------------------------//
bool Element_Definition::invariant_satisfied() const
{
    bool ldum = ( name.empty() == false );
    
    if( type == NODE )
    {
	ldum = ldum && (dimension         == 0);
	ldum = ldum && (number_of_nodes   == 1);
	ldum = ldum && (number_of_sides   == 0);	
	ldum = ldum && (elem_defs.size()  == 0);
    }
    else
    {
	ldum = ldum && (dimension > 0); 
	ldum = ldum && (dimension < 4);
	ldum = ldum && (number_of_nodes > dimension);
	ldum = ldum && (number_of_sides <= number_of_nodes);
	ldum = ldum && (number_of_sides > dimension); 
	ldum = ldum && (elem_defs.size() > 0);
    }
    
    ldum = ldum && (side_type.size()  == number_of_sides);
    ldum = ldum && (side_nodes.size() == number_of_sides);
    ldum = ldum && (node_loc.size()   == number_of_nodes);
    
    for( size_t i=0; i < elem_defs.size(); i++ )
	ldum = ldum && (elem_defs[i].dimension == dimension-1);
    
    for( size_t i=0; i < side_nodes.size(); i++ )
    {
	ldum = ldum && (side_nodes[i].size() > 0);
	ldum = ldum && (side_nodes[i].size() == 
			elem_defs[ side_type[i] ].number_of_nodes);
	for( size_t j=0; j < side_nodes[i].size(); j++ )
	{
	    // ldum = ldum && (side_nodes[i][j] >= 0);
	    ldum = ldum && (side_nodes[i][j] < number_of_nodes);
	    ldum = ldum && (node_loc[ side_nodes[i][j] ] ==
			    elem_defs[ side_type[i] ].node_loc[j]);
	}
    }

    return ldum;
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_node()
{
    name = "NODE";
    dimension=0;
    number_of_sides=0;
    number_of_nodes=1;
    node_loc.push_back(CORNER);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_bar()
{
    std::vector<size_t> tmp;
    dimension=1;
    number_of_sides=2;
    tmp.clear();
    tmp.push_back(0);
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    side_nodes.push_back(tmp);
    for( size_t i=0; i < 2; i++ )
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
	Insist(false,"#2 Unrecognized Element-Type Flag");
    }
    elem_defs.push_back( Element_Definition(NODE) );
    for( size_t i = 0; i < number_of_sides; i++ )
	side_type.push_back(0);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_tri()
{
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 3; i++ )
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
	for (size_t i=0; i < 3; i++)
	    node_loc.push_back(EDGE);
	break;
    default :
	Insist(false,"#3 Unrecognized Element-Type Flag");
    } 

    for( size_t i = 0; i < number_of_sides; i++ )
	side_type.push_back(0);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_quad()
{
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 4; i++ )
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
	for (size_t i=0; i < 4; i++)
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
	    Insist(false,"#4 Unrecognized Element-Type Flag");
	}
	break;
    default :
	Insist(false,"#5 Unrecognized Element-Type Flag");
    }

    for( size_t i = 0; i < number_of_sides; i++ )
	side_type.push_back(0);
}

//---------------------------------------------------------------------------//

/*
void Element_Definition::construct_pentagon()
{
    return;
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 4; i++ )
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
        for (size_t i=0; i < 4; i++)
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
            Insist(false,"#4 Unrecognized Element-Type Flag");
        }
        break;
    default :
        Insist(false,"#5 Unrecognized Element-Type Flag");
    }

    for( size_t i = 0; i < number_of_sides; i++ )
        side_type.push_back(0);
}
*/

//---------------------------------------------------------------------------//

void Element_Definition::construct_pentagon()
{
    std::vector<size_t> tmp;
    dimension=2;
    number_of_sides=5;
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
    tmp[1] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 4;
    tmp[1] = 0;
    side_nodes.push_back(tmp);
    for( size_t i=0; i < 4; i++ )
	node_loc.push_back(CORNER);

    switch ( type )
    {
    case PENTAGON_5 :
	name = "PENTAGON_5";
	number_of_nodes=5;
	elem_defs.push_back(Element_Definition(BAR_2));
	break;
    default :
	Insist(false,"#5 Unrecognized Element-Type Flag");
    }

    for( size_t i = 0; i < number_of_sides; i++ )
	side_type.push_back(0);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_tetra()
{
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 4; i++ )
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
	
	for (size_t i=0; i < 6; i++)
	    node_loc.push_back(EDGE);
	break;
    default :
	Insist(false,"#6 Unrecognized Element-Type Flag");
    }
    for (size_t i = 0; i < number_of_sides; i++)
	side_type.push_back(0);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_pyra()
{
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 5; i++ )
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
	for (size_t i=0; i < 8; i++)
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
	Insist(false,"#7 Unrecognized Element-Type Flag");
    }

    side_type.push_back(0);
    for( size_t i = 1; i < number_of_sides; i++ )
	side_type.push_back(1);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_penta()
{
    std::vector<size_t> tmp;
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
    for( size_t i=0; i < 6; i++ )
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
	for (size_t i=0; i < 9; i++)
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
	    for (size_t i=0; i < 3; i++)
		node_loc.push_back(FACE);
	    number_of_nodes=18;
	    elem_defs.push_back(Element_Definition(QUAD_9));
	    elem_defs.push_back(Element_Definition(TRI_6));
	    side_nodes[0].push_back(15);
	    side_nodes[1].push_back(16);
	    side_nodes[2].push_back(17);
	    break;
	default :
	    Insist(false,"#8 Unrecognized Element-Type Flag");
	}
	break;
    default :
	Insist(false,"#9 Unrecognized Element-Type Flag");
    }

    for( size_t i = 0; i < 3; i++ )
	side_type.push_back(0);
    for( size_t i = 3; i < number_of_sides; i++ )
	side_type.push_back(1);
}

//---------------------------------------------------------------------------//

void Element_Definition::construct_hexa()
{
    std::vector<size_t> tmp;
    dimension=3;
    number_of_sides=6;
    tmp.clear();
    tmp.push_back(0);
    tmp.push_back(3); 
    tmp.push_back(2); 
    tmp.push_back(1);   
    side_nodes.push_back(tmp);
    tmp[0] = 0;
    tmp[1] = 4;
    tmp[2] = 7;
    tmp[3] = 3;
    side_nodes.push_back(tmp);
    tmp[0] = 2;
    tmp[1] = 3;
    tmp[2] = 7;
    tmp[3] = 6;
    side_nodes.push_back(tmp);
    tmp[0] = 1;
    tmp[1] = 2;
    tmp[2] = 6;
    tmp[3] = 5;
    side_nodes.push_back(tmp);
    tmp[0] = 0;
    tmp[1] = 1;
    tmp[2] = 5;
    tmp[3] = 4;
    side_nodes.push_back(tmp);
    tmp[0] = 4;
    tmp[1] = 5;
    tmp[2] = 6;
    tmp[3] = 7;
    side_nodes.push_back(tmp);
    for( size_t i=0; i < 8; i++ )
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
	for (size_t i=0; i < 12; i++)
	    node_loc.push_back(EDGE);
	
	side_nodes[0].push_back(11);
	side_nodes[0].push_back(10);
	side_nodes[0].push_back(9);
	side_nodes[0].push_back(8);
	
	side_nodes[1].push_back(12);
	side_nodes[1].push_back(19);
	side_nodes[1].push_back(15);
	side_nodes[1].push_back(11);
	
	side_nodes[2].push_back(10);
	side_nodes[2].push_back(15);
	side_nodes[2].push_back(18);
	side_nodes[2].push_back(14);
	
	side_nodes[3].push_back(9);
	side_nodes[3].push_back(14);
	side_nodes[3].push_back(17);
	side_nodes[3].push_back(13);
	
	side_nodes[4].push_back(8);
	side_nodes[4].push_back(13);
	side_nodes[4].push_back(16);
	side_nodes[4].push_back(12);
	
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
	    for (size_t i=0; i < 6; i++)
		node_loc.push_back(FACE);
	    node_loc.push_back(CELL);
	    number_of_nodes=27;
	    elem_defs.push_back(Element_Definition(QUAD_9));
	    for (size_t i=0; i < number_of_sides; i++)
		side_nodes[i].push_back(i+20);
	    break;
	default :
	    Insist(false,"#10 Unrecognized Element-Type Flag");
	}
	break;
    default :
	Insist(false,"#11 Unrecognized Element-Type Flag");
    }

    for( size_t i = 0; i < number_of_sides; i++ )
	side_type.push_back(0);	
}

//---------------------------------------------------------------------------//

std::ostream& Element_Definition::print( std::ostream & os_out ) const
{
    os_out << "Element Type   : " << get_type() << std::endl;
    os_out << "Element Name   : " << get_name() << std::endl;
    os_out << "Number of Nodes: " << get_number_of_nodes() <<
	std::endl;
    os_out << "Dimension      : " << get_dimension() << std::endl;
    os_out << "Number of Sides: " << get_number_of_sides() << std::endl;
    os_out << "Node Locations : ";
    for( size_t j=0; j<get_number_of_nodes(); j++ )
	os_out << get_node_location(j) << " ";
    os_out << std::endl;
    if( get_number_of_sides() != 0 ) 
    {
	os_out << "Side Types     : ";
	for( size_t j=0; j<get_number_of_sides(); j++ )
	    os_out << get_side_type(j).get_name() << " ";
	os_out << std::endl;
	
	std::vector< size_t > tmp;
	os_out << "Side Nodes     : " << std::endl;
	for( size_t j=0; j<get_number_of_sides(); j++ )
	{
	    tmp = get_side_nodes(j);
	    os_out << "  " << "side# " << j << " -    ";
	    for( size_t k=0; k<tmp.size(); k++ )
		os_out << tmp[k] << " ";
	    os_out << std::endl;
	}
    }
    {
        std::vector<unsigned> num_face_nodes = get_number_of_face_nodes();
        std::vector<std::vector<unsigned> > face_nodes = get_face_nodes();
        os_out <<"Face Nodes: " << std::endl;
        for (size_t j=0; j<num_face_nodes.size(); ++j)
        {
            os_out << "  Face " << j << ": " << num_face_nodes[j]
                   << " nodes : ";
            for (size_t k=0; k<num_face_nodes[j]; ++k)
            {
                os_out << face_nodes[j][k] << " ";
            }
            os_out << std::endl;
        }
    }
    os_out << std::endl;
    return os_out;
}

} // end namespace rtt_mesh_element

//---------------------------------------------------------------------------//
//                       end of Element_Definition.cc
//---------------------------------------------------------------------------//
