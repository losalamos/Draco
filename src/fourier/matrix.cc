//----------------------------------*-C++-*----------------------------------//
// matrix.cc
// John Gulick
// Fri Aug  6 10:55:21 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matrix.t.hh"


//**************************
//Constructor function body.
//**************************
//Constructor
Matrix::Matrix(int mat_order)
{
this->mat_order        = mat_order;
T_pos = new complex<double> [2*mat_order];
T_neg = new complex<double> [2*mat_order];
P     = new complex<double> [2*mat_order];
H     = new complex<double> [2*mat_order];
}


//********************************
//Matrix building member function.
//********************************
//SCB matrix building member function.
lapack_matrix<complex<double> >::type Matrix::SCB_matrix(double sigma_t,
						      double sigma_s,
						      int quad_order, 
						      double lambda, 
						      double delta_x)  

{
//Setup quadrature.
Gauss quadrature;
quadrature.gauleg(quad_order);
mu     = quadrature.get_mu();
weight = quadrature.get_weight();
gamma  = quadrature.get_gamma();

typedef complex<double> Type;

lapack_matrix<Type>::type T_positive(mat_order,mat_order);
lapack_matrix<Type>::type T_negative(mat_order,mat_order);
lapack_matrix<Type>::type T_total(mat_order,mat_order);
for(int j=1; j<=quad_order/2; j++)
    { 
    //Create T positive matrix.
    T_pos[0] = Type (((mu[j]/2.) + (sigma_t*delta_x/2.0)),(0.));
    T_pos[2] = Type (((mu[j]/2.) - (mu[j]*cos(lambda))),(mu[j]*sin(lambda)));
    T_pos[1] = Type ((-mu[j]/2.),(0.));
    T_pos[3] = Type ((mu[j] - (mu[j]/2.) + (sigma_t*delta_x/2.)),(0.));

    //Create T negative matrix.
    T_neg[0] = Type ((mu[j] - (mu[j]/2.) + (sigma_t*delta_x/2.)),(0.));
    T_neg[2] = Type ((-mu[j]/2.),(0.));
    T_neg[1] = Type (((mu[j]/2.) - (mu[j]*cos(lambda))),(-mu[j]*sin(lambda)));
    T_neg[3] = Type (((mu[j]/2.) + (sigma_t*delta_x/2.0)),(0.));


    //Create P (source) matrix.
    P[0] = Type ((sigma_s*delta_x/4.),(0.));
    P[2] = Type ((0.),(0.));
    P[1] = Type ((0.),(0.));
    P[3] = Type ((sigma_s*delta_x/4.),(0.));

    //Define T_pos, T_neg, P matrices as lapack_matrix type.
    lapack_matrix<Type,external>::type T_pos_mat(T_pos,mat_order,mat_order);
    lapack_matrix<Type,external>::type T_neg_mat(T_neg,mat_order,mat_order);
    lapack_matrix<Type,external>::type P_mat(P,mat_order,mat_order);

    //Take the inverse of the T_pos and T_neg matrices.
    lapack_matrix<Type>::type B_pos(mat_order,mat_order);
    lapack_matrix<Type>::type B_neg(mat_order,mat_order);
    mtl::set(B_pos,0);
    mtl::set(B_neg,0);
    set_diagonal(B_pos,1);
    set_diagonal(B_neg,1);
    dense1D<int> ipivot_pos(mat_order,0);
    dense1D<int> ipivot_neg(mat_order,0);
    int info_pos;
    int info_neg;
	
    //T_positive inversion.
    info_pos = getrf(T_pos_mat,ipivot_pos);
    if(info_pos == 0)
       {
       info_pos = getrs('N',T_pos_mat,ipivot_pos,B_pos);
       if(info_pos != 0)
	  {
	  cerr << "Inversion failed at getrs(pos)." << info_pos <<endl;
	  }
       }
    else
       {
       cerr << "Inversion failed at getrf(pos)." << info_pos <<endl;
       }
	
    //T_negative inversion.
    info_neg = getrf(T_neg_mat,ipivot_neg);
    if(info_neg == 0)
       {
       info_neg = getrs('N',T_neg_mat,ipivot_neg,B_neg);
       if(info_neg != 0)
	  {
	  cerr << "Inversion failed at getrs(neg)." << info_neg <<endl;
	  }
       }
    else
       {
       cerr << "Inversion failed at getrf(neg)." << info_neg <<endl;
       }       

    //Multiply inverse T matrix by P to get T_(pos/neg)_tot.
    lapack_matrix<Type>::type T_pos_tot(mat_order,mat_order);
    lapack_matrix<Type>::type T_neg_tot(mat_order,mat_order);
    mtl::set(T_pos_tot,0.);
    mtl::set(T_neg_tot,0.);
    mult(B_pos,P_mat,T_pos_tot);
    mult(B_neg,P_mat,T_neg_tot);

    //Summation of T_positive and T_negative.
    scale(T_pos_tot,weight[j]);
    scale(T_neg_tot,weight[j]);
    add(T_pos_tot,T_positive);
    add(T_neg_tot,T_negative);
    }

//Summation to create T_total matrix.
add(T_positive,T_total);
add(T_negative,T_total);

//Return T_total.  The A matrix in the eigensystem A*x=lam*x.
return T_total;
}


//---------------------------------------------------------------------------//
//                              end of matrix.cc
//---------------------------------------------------------------------------//
