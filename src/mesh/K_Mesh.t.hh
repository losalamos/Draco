//----------------------------------*-C++-*----------------------------------//
// K_Mesh.cc
// Geoffrey Furnish
// Wed Jun 11 00:14:08 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "K_Mesh.hh"

template<class Geometry>
K_Mesh<Geometry>::K_Mesh( const K_Mesh_DB& kdb )
    : Geometry(),
      K_Mesh_DB( kdb ),
      C( ncx ), F( ncx+1 ),
      clayout( C ), flayout( F ),
      xc( this ), Vc( this ),
      xf( this ), Af( this )
{
    cout << "Constructing a K_Mesh" << endl;
    cout << "ncx=" << ncx << endl;

// Now need to set up the elements.

    xf[F] = F*(xmax-xmin)/(ncx-1);
}

//---------------------------------------------------------------------------//
// Define the nested class stuff.

template<class Geometry>
K_Mesh<Geometry>::ccsf::ccsf( K_Mesh<Geometry> *pm )
    : Field<double, 1>( pm->clayout )
{
    cout << "ccsf, bare * ctor" << endl;
}

template<class Geometry>
K_Mesh<Geometry>::ccsf::ccsf( SP< K_Mesh<Geometry> >& m_ )
    : Field<double, 1>( m_->clayout ),
      m(m_)
{
    cout << "ccsf, sp ctor" << endl;
}

template<class Geometry>
K_Mesh<Geometry>::fcsf::fcsf( K_Mesh<Geometry> *pm )
    : Field<double, 1>( pm->flayout )
{
    cout << "fcsf, bare * ctor" << endl;
}

template<class Geometry>
K_Mesh<Geometry>::fcsf::fcsf( SP< K_Mesh<Geometry> >& m_ )
    : Field<double, 1>( m_->flayout ),
      m(m_)
{
    cout << "fcsf, sp ctor" << endl;
}

//---------------------------------------------------------------------------//
//                              end of K_Mesh.cc
//---------------------------------------------------------------------------//
