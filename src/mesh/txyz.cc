//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Geoffrey M. Furnish
// Wed May 13 09:53:44 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
using namespace std;

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    cout << "t1: passed\n";

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    dsxx::SP<Mesh_XYZ> spm = new Mesh_XYZ( mdb );

    cout << "t2: passed" << endl;

    {
        Mesh_XYZ* m1 = new Mesh_XYZ( mdb );
        Mesh_XYZ* m2 = new Mesh_XYZ( mdb );
        if (*m1 == *m2)
            cout << "error: same mesh" << endl;
        m1 = m2;
        if (!(*m1 == *m2))
            cout << "error: different meshes" << endl;
    }

    {
        Mesh_XYZ::ccsf x( spm ), y( spm ), z( spm );
        x = 1.;
        y = 2.;
        z = x + y;
        Mesh_XYZ m = x.get_Mesh();
    }

    {
        Mesh_XYZ::ccif xi( spm ), yi( spm ), zi( spm );
        xi = 1;
        yi = 2;
        zi = xi + yi;
        Mesh_XYZ m = xi.get_Mesh();
    }

    {
        Mesh_XYZ::fcdsf xf( spm ), yf( spm ), zf( spm );
        xf = 1.;
        yf = xf;
        zf = xf + yf;
        zf = xf - yf;
        Mesh_XYZ m = xf.get_Mesh();
    }

    {
        Mesh_XYZ::ncsf xn( spm ), yn( spm ), zn( spm );
        xn = 1.;
        yn = 2.;
        zn = xn + yn;
        Mesh_XYZ m = xn.get_Mesh();
    }

    {
        Mesh_XYZ::vcsf xv( spm ), yv( spm ), zv( spm );
        xv = 1.;
        yv = 2.;
        zv = xv + yv;
        Mesh_XYZ m = xv.get_Mesh();
    }

    {
        Mesh_XYZ::vec3 xvec, yvec, zvec;
        double dotprod;
        xvec = 1.;
        yvec = 2.;
        zvec = xvec + yvec;
        dotprod = Mesh_XYZ::vec3::dot(xvec, yvec);
	if (dotprod < 5.999 || dotprod > 6.001)
            cout << "Error in dot product" << endl;
    }

    {
        Mesh_XYZ::gccsf xgc( spm );
        Mesh_XYZ::ccsf x( spm );
        x = 1.;
        xgc = x;
        Mesh_XYZ m = xgc.get_Mesh();
    }

    {
        Mesh_XYZ::gfcdsf xgf( spm );
        Mesh_XYZ::fcdsf xf( spm );
        xf = 1.;
        xgf = xf;
        Mesh_XYZ m = xgf.get_Mesh();
    }

    {
        Mesh_XYZ::gvcsf xgv( spm );
        Mesh_XYZ::vcsf xv( spm );
        xv = 1.;
        xgv = xv;
        Mesh_XYZ m = xgv.get_Mesh();
    }

    {
      try
      {
        Mesh_XYZ::bssf xb( spm ), yb( spm ), zb( spm );
        xb = 1.;
        yb = 2.;
        zb = xb + yb;
        Mesh_XYZ m = xb.get_Mesh();
        for (Mesh_XYZ::bssf::iterator iter = xb.begin();
             iter != xb.end(); ++iter)
            *iter = 6.;
        for (Mesh_XYZ::bssf::const_iterator iter = zb.begin();
             iter != zb.end(); ++iter)
	    if (*iter < 2.999 || *iter > 3.001)
                cout << "Error in boundary treatment" << endl;

        const Mesh_XYZ::bssf cxb( spm );
        Mesh_XYZ::bssf::iterator iter = xb.begin();
        for (Mesh_XYZ::bssf::const_iterator citer = cxb.begin();
             citer != cxb.end(); ++citer, ++iter)
            *iter = *citer;
      }
      catch (const dsxx::assertion &ass)
      {
        std::cerr << "Caught dsxx::assertion bssf: '"
	          << ass.what() << "'." << endl;
      }
    }

    {
        Mesh_XYZ::ccsf oneCC( spm );
        oneCC = 1.0;
        //dump( oneCC, "oneCC, before" );
        Mesh_XYZ::fcdsf twoFC( spm );
        twoFC = 0.0;
        //dump( twoFC, "twoFC, before" );
        Mesh_XYZ::scatter ( twoFC, oneCC, Mesh_XYZ::OpAddAssign() );
        //dump( oneCC, "oneCC, after" );
        //dump( twoFC, "twoFC, after" );
    }

    {
        Mesh_XYZ::ccsf threeCC( spm );
        threeCC = 3.0;
        //dump( threeCC, "threeCC, before" );
        Mesh_XYZ::fcdsf nineFC( spm );
        nineFC = 1.0;
        //dump( nineFC, "nineFC, before" );
        Mesh_XYZ::scatter ( nineFC, threeCC, Mesh_XYZ::OpMultAssign() );
        //dump( threeCC, "threeCC, after" );
        //dump( nineFC, "nineFC, after" );
    }

    {
        Mesh_XYZ::fcdsf oneFC( spm );
        oneFC = 1.0;
        Mesh_XYZ::ccsf sixCC( spm );
        sixCC = 0.0;
        Mesh_XYZ::scatter ( sixCC, oneFC, Mesh_XYZ::OpAddAssign() );
        //dump( sixCC, "sixCC, after" );
    }

    {
        Mesh_XYZ::fcdsf threeFC( spm );
        threeFC = 3.0;
        Mesh_XYZ::ccsf CC729( spm );
        CC729 = 1.0;
        Mesh_XYZ::scatter ( CC729, threeFC, Mesh_XYZ::OpMultAssign() );
        //dump( CC729, "CC729, after" );
    }

    {
        Mesh_XYZ::fcdsf oneFC( spm );
        oneFC = 0.0;
        Mesh_XYZ::ccsf oneCC( spm );
        oneCC = 1.0;
        Mesh_XYZ::gather ( oneFC, oneCC, Mesh_XYZ::OpAddAssign() );
        //dump( oneFC, "oneFC, after" );
    }

    {
        Mesh_XYZ::vcsf oneVC( spm );
        oneVC = 0.0;
        Mesh_XYZ::ncsf oneNC( spm );
        oneNC = 1.0;
        Mesh_XYZ::gather ( oneVC, oneNC, Mesh_XYZ::OpAddAssign() );
        //dump( oneNC, "oneNC, after" );
        //dump( oneVC, "oneVC, after" );
    }

    {
        Mesh_XYZ::ncsf eightNC( spm );
        eightNC = 0.0;
        Mesh_XYZ::vcsf oneVC( spm );
        oneVC = 1.0;
        Mesh_XYZ::scatter ( eightNC, oneVC, Mesh_XYZ::OpAddAssign() );
        //dump( eightNC, "eightNC, after" );
    }

    {
        Mesh_XYZ::vcsf oneVC( spm );
        oneVC = 1.0;
        Mesh_XYZ::fcdsf fourFC( spm );
        fourFC = 0.0;
        Mesh_XYZ::scatter ( fourFC, oneVC, Mesh_XYZ::OpAddAssign() );
        //dump( fourFC, "fourFC, after" );
    }

    {
        Mesh_XYZ::vcsf threeVC( spm );
        threeVC = 3.0;
        Mesh_XYZ::fcdsf FC81( spm );
        FC81 = 1.0;
        Mesh_XYZ::scatter ( FC81, threeVC, Mesh_XYZ::OpMultAssign() );
        //dump( FC81, "FC81, after" );
    }

    {
        Mesh_XYZ::vcvf oneVCV( spm );
        Mesh_XYZ::vec3 onevec;
        onevec = 1.0;
        oneVCV = onevec;
        Mesh_XYZ::fcdvf fourFCV( spm );
        Mesh_XYZ::vec3 zerovec;
        zerovec = 0.0;
        fourFCV = zerovec;
        Mesh_XYZ::scatter( fourFCV, oneVCV, Mesh_XYZ::OpAddAssign() );
    }

    {
        Mesh_XYZ::fcdsf swapFC( spm );
        swapFC = 0.;
        Mesh_XYZ::fcdsf oneFC( spm );
        oneFC = 1.;
        Mesh_XYZ::swap( swapFC, oneFC );
        //dump( swapFC, "swapFC, after" );
    }

    {
        Mesh_XYZ::ccsf oneCC( spm );
        oneCC = 1.;
        double total;
        total = Mesh_XYZ::sum(oneCC);
        if (total < 63.999 || total > 64.001)
            cout << "Error in global sum" << endl;
    }

    {
        Mesh_XYZ::ccif oneCC( spm );
        oneCC = 1;
        int total;
        total = Mesh_XYZ::sum(oneCC);
        if (total != 64)
            cout << "Error in global sum" << endl;
    }

    {
        Mesh_XYZ::fcdsf oneFC( spm );
        oneFC = 1.;
        double total;
        total = Mesh_XYZ::sum(oneFC);
        if (total < 383.999 || total > 384.001)
            cout << "Error in global sum" << endl;
    }

    {
        Mesh_XYZ::bssf oneBS( spm );
        oneBS = 0.;
        Mesh_XYZ::fcdsf oneFC( spm );
        oneFC = 1.;
        Mesh_XYZ::gather( oneBS, oneFC, Mesh_XYZ::OpAssign() );
        oneFC = 0.;
        Mesh_XYZ::gather( oneFC, oneBS, Mesh_XYZ::OpAssign() );
        //dump( oneFC, "oneFC, after" );
    }

    {
        Mesh_XYZ::fcdvf face_normals( spm );
        face_normals = spm->get_fn();

        Mesh_XYZ::fcdsf face_areas( spm );
        spm->get_face_areas(face_areas);

        Mesh_XYZ::fcdsf face_lengths( spm );
        spm->get_face_lengths(face_lengths);

        Mesh_XYZ::fcdsf face_locs( spm );
        spm->get_xloc(face_locs);
        spm->get_yloc(face_locs);
        spm->get_zloc(face_locs);
    }

    {
        Mesh_XYZ::ccsf cell_locs( spm );
        spm->get_xloc(cell_locs);
        spm->get_yloc(cell_locs);
        spm->get_zloc(cell_locs);

        Mesh_XYZ::ccsf cell_deltas( spm );
        spm->get_dx(cell_deltas);
        spm->get_dy(cell_deltas);
        spm->get_dz(cell_deltas);

      try
      {
        Mesh_XYZ::ccsf cell_volumes( spm );
        spm->get_cell_volumes(cell_volumes);
      }
      catch (const dsxx::assertion &ass)
      {
        std::cerr << "Caught dsxx::assertion cell_volumes: '"
	          << ass.what() << "'." << endl;
      }
    }

    cout << "Finalizing" << endl;

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
