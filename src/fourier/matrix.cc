//----------------------------------*-C++-*----------------------------------//
// matrix.cc
// John Gulick
// Fri Aug  6 10:55:21 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matrix.hh"


//**************************
//Constructor function body.
//**************************
//Matrix class constructor.
Matrix::Matrix(int mat_order)
{
    this->mat_order = mat_order;
}


//********************************
//Matrix building member function.
//********************************
//Simple Corner Balance (SCB) matrix building member function.
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

    T_pos = new complex<double> [2*mat_order];
    T_neg = new complex<double> [2*mat_order];
    P     = new complex<double> [2*mat_order];

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


//Simple Corner Balance (SCB) accelerated with SCB derived "Modified 4-Step"
//DSA  matrix building member function.
lapack_matrix<complex<double> >::type Matrix::SCB_M4S_matrix(double sigma_t,
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


    T_pos = new complex<double> [2*mat_order];
    T_neg = new complex<double> [2*mat_order];
    P     = new complex<double> [2*mat_order];
    H     = new complex<double> [2*mat_order];

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

    //Create the SCB derived DSA matrix.
    double c         = sigma_s/sigma_t;
    double one       = 1./(3.*sigma_t*sigma_t*delta_x*delta_x);
    double two       = gamma/(sigma_t*delta_x);
    double three     = 1-c;

    H[0] = Type ((one + two + three - one*cos(lambda)),(one*sin(lambda)));
    H[1] = Type ((-one + (one-two)*cos(lambda)),(-(two-one)*sin(lambda)));
    H[2] = Type ((-one + (one-two)*cos(lambda)),(-(one-two)*sin(lambda)));
    H[3] = Type ((one + two + three - one*cos(lambda)),(-one*sin(lambda)));
    lapack_matrix<Type,external>::type H_dsa(H,mat_order,mat_order);

    //H_dsa inversion.
    lapack_matrix<Type>::type B_dsa(mat_order,mat_order);
    mtl::set(B_dsa,0);
    set_diagonal(B_dsa,1);
    dense1D<int> ipivot_dsa(mat_order,0);
    int info_dsa;
    info_dsa = getrf(H_dsa,ipivot_dsa);
    if(info_dsa == 0)
    {
	info_dsa = getrs('N',H_dsa,ipivot_dsa,B_dsa);
	if(info_dsa != 0)
	{
	    cerr << "Inversion failed at getrs(dsa)." << info_dsa <<endl;
	}
    }
    else
    {
	cerr << "Inversion failed at getrf(dsa)." << info_dsa <<endl;
    }

    //Multiply the inverted matrix B_dsa by c.
    scale(B_dsa,c);
    lapack_matrix<Type>::type ID(mat_order,mat_order);
    lapack_matrix<Type>::type T_total_tmp(mat_order,mat_order);
    copy(T_total,T_total_tmp);
    set_diagonal(ID,-1);
    add(ID,T_total);
    lapack_matrix<Type>::type T_total_mat(mat_order,mat_order);
    mtl::set(T_total_mat,0);
    mult(T_total,B_dsa,T_total_mat);
    add(T_total_mat,T_total_tmp);
    T_total = T_total_tmp;

    //Return T_total.  The A matrix in the eigensystem A*x=lam*x.
    return T_total;
}


//Upstream Corner Balance (UCB) matrix building member function.
lapack_matrix<complex<double> >::type Matrix::UCB_matrix(double sigma_t,
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

    T_pos = new complex<double> [2*mat_order];
    T_neg = new complex<double> [2*mat_order];
    P     = new complex<double> [2*mat_order];

    lapack_matrix<Type>::type T_positive(mat_order,mat_order);
    lapack_matrix<Type>::type T_negative(mat_order,mat_order);
    lapack_matrix<Type>::type T_total(mat_order,mat_order);
    for(int j=1; j<=quad_order/2; j++)
    { 
	//Determine UCB coefficients.
	double alpha_not = 0.455;
	double theta     = 1.;
	double c         = sigma_s/sigma_t;
	double tau       = sigma_t* delta_x/(2.*mu[j]);
	double alpha     = (3. + (4.*tau) + (alpha_not*4.*tau*tau))/
	    (2. + (2.*tau) + (4.*tau*tau));
	double beta      = alpha/tau;

	//Create T positive matrix.
        T_pos[0] = Type (((mu[j]/delta_x)+((mu[j]*beta)/delta_x)+(1./2.)),(0.));
	T_pos[2] = Type (((mu[j]/delta_x)*(cos(-lambda*sigma_t)*(-1.-beta))),
			 ((mu[j]/delta_x)*(sin(-lambda*sigma_t))*(-1.-beta)));
	T_pos[1] = Type (((mu[j]/delta_x)*(((1-theta)/2.)-1.-beta)),(0.)); 
	T_pos[3] = Type (((mu[j]*(1.+theta)/(2.*delta_x))+
			  (mu[j]*beta/delta_x)*(cos(-lambda*sigma_t))+(1./2.)),
			 ((mu[j]*beta/delta_x)*sin(-lambda*sigma_t)));

	//Create T negative matrix.
	T_neg[0] = Type (((mu[j]*(1.+theta)/(2.*delta_x))+
			  (mu[j]*beta/delta_x)*(cos(lambda*sigma_t)) + (1./2.)),
			 ((mu[j]*beta/delta_x)*sin(lambda*sigma_t)));
	T_neg[2] = Type (((mu[j]/delta_x)*(((1.-theta)/2.)-1.-beta)),(0.));
	T_neg[1] = Type (((mu[j]/delta_x)*(cos(lambda*sigma_t)*(-1.-beta))),
			 ((mu[j]/delta_x)*(sin(lambda*sigma_t))*(-1.-beta)));
        T_neg[3] = Type (((mu[j]/delta_x)+((mu[j]*beta)/delta_x)+(1./2.)),(0.0));


	//Create P (source) matrix.
	P[0] = Type ((c/4)*(1+(mu[j]/delta_x)),(0.0));
	P[2] = Type ((c/4)*( -(mu[j]/delta_x)),(0.0));
	P[1] = Type ((c/4)*( -(mu[j]/delta_x)),(0.0));
	P[3] = Type ((c/4)*(1+(mu[j]/delta_x)),(0.0));

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

//Upstream Corner Balance (UCB) accelerated with SCB derived "Modified 4-Step"
//DSA  matrix building member function.
lapack_matrix<complex<double> >::type Matrix::UCB_M4S_matrix(double sigma_t,
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


    T_pos = new complex<double> [2*mat_order];
    T_neg = new complex<double> [2*mat_order];
    P     = new complex<double> [2*mat_order];
    H     = new complex<double> [2*mat_order];

    lapack_matrix<Type>::type T_positive(mat_order,mat_order);
    lapack_matrix<Type>::type T_negative(mat_order,mat_order);
    lapack_matrix<Type>::type T_total(mat_order,mat_order);
    for(int j=1; j<=quad_order/2; j++)
    { 
	//Determine UCB coefficients.
	double alpha_not = 0.455;
	double theta     = 1.;
	double c         = sigma_s/sigma_t;
	double tau       = sigma_t* delta_x/(2.*mu[j]);
	double alpha     = (3. + (4.*tau) + (alpha_not*4.*tau*tau))/
	    (2. + (2.*tau) + (4.*tau*tau));
	double beta      = alpha/tau;

	//Create T positive matrix.
        T_pos[0] = Type (((mu[j]/delta_x)+((mu[j]*beta)/delta_x)+(1./2.)),(0.));
	T_pos[2] = Type (((mu[j]/delta_x)*(cos(-lambda*sigma_t)*(-1.-beta))),
			 ((mu[j]/delta_x)*(sin(-lambda*sigma_t))*(-1.-beta)));
	T_pos[1] = Type (((mu[j]/delta_x)*(((1-theta)/2.)-1.-beta)),(0.)); 
	T_pos[3] = Type (((mu[j]*(1.+theta)/(2.*delta_x))+
			  (mu[j]*beta/delta_x)*(cos(-lambda*sigma_t))+(1./2.)),
			 ((mu[j]*beta/delta_x)*sin(-lambda*sigma_t)));

	//Create T negative matrix.
	T_neg[0] = Type (((mu[j]*(1.+theta)/(2.*delta_x))+
			  (mu[j]*beta/delta_x)*(cos(lambda*sigma_t)) + (1./2.)),
			 ((mu[j]*beta/delta_x)*sin(lambda*sigma_t)));
	T_neg[2] = Type (((mu[j]/delta_x)*(((1.-theta)/2.)-1.-beta)),(0.));
	T_neg[1] = Type (((mu[j]/delta_x)*(cos(lambda*sigma_t)*(-1.-beta))),
			 ((mu[j]/delta_x)*(sin(lambda*sigma_t))*(-1.-beta)));
        T_neg[3] = Type (((mu[j]/delta_x)+((mu[j]*beta)/delta_x)+(1./2.)),(0.0));


	//Create P (source) matrix.
	P[0] = Type ((c/4)*(1+(mu[j]/delta_x)),(0.0));
	P[2] = Type ((c/4)*( -(mu[j]/delta_x)),(0.0));
	P[1] = Type ((c/4)*( -(mu[j]/delta_x)),(0.0));
	P[3] = Type ((c/4)*(1+(mu[j]/delta_x)),(0.0));

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

    //Create the SCB derived DSA matrix.
    double c         = sigma_s/sigma_t;
    double one       = 1./(3.*sigma_t*sigma_t*delta_x*delta_x);
    double two       = gamma/(sigma_t*delta_x);
    double three     = 1-c;

    H[0] = Type ((one + two + three - one*cos(lambda)),(one*sin(lambda)));
    H[1] = Type ((-one + (one-two)*cos(lambda)),(-(two-one)*sin(lambda)));
    H[2] = Type ((-one + (one-two)*cos(lambda)),(-(one-two)*sin(lambda)));
    H[3] = Type ((one + two + three - one*cos(lambda)),(-one*sin(lambda)));
    lapack_matrix<Type,external>::type H_dsa(H,mat_order,mat_order);

    //H_dsa inversion.
    lapack_matrix<Type>::type B_dsa(mat_order,mat_order);
    mtl::set(B_dsa,0);
    set_diagonal(B_dsa,1);
    dense1D<int> ipivot_dsa(mat_order,0);
    int info_dsa;
    info_dsa = getrf(H_dsa,ipivot_dsa);
    if(info_dsa == 0)
    {
	info_dsa = getrs('N',H_dsa,ipivot_dsa,B_dsa);
	if(info_dsa != 0)
	{
	    cerr << "Inversion failed at getrs(dsa)." << info_dsa <<endl;
	}
    }
    else
    {
	cerr << "Inversion failed at getrf(dsa)." << info_dsa <<endl;
    }

    //Multiply the inverted matrix B_dsa by c.
    scale(B_dsa,c);
    lapack_matrix<Type>::type ID(mat_order,mat_order);
    lapack_matrix<Type>::type T_total_tmp(mat_order,mat_order);
    copy(T_total,T_total_tmp);
    set_diagonal(ID,-1);
    add(ID,T_total);
    lapack_matrix<Type>::type T_total_mat(mat_order,mat_order);
    mtl::set(T_total_mat,0);
    mult(T_total,B_dsa,T_total_mat);
    add(T_total_mat,T_total_tmp);
    T_total = T_total_tmp;

    //Return T_total.  The A matrix in the eigensystem A*x=lam*x.
    return T_total;
}


//Linear Characteristics Method matrix building routine.
lapack_matrix<complex<double> >::type Matrix::LIN_CHAR_matrix(double sigma_t,
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

    T_pos = new complex<double> [2*mat_order];
    T_neg = new complex<double> [2*mat_order];
    P     = new complex<double> [2*mat_order];

    lapack_matrix<Type>::type T_positive(mat_order,mat_order);
    lapack_matrix<Type>::type T_negative(mat_order,mat_order);
    lapack_matrix<Type>::type T_total(mat_order,mat_order);
    for(int j=1; j<=quad_order/2; j++)
    { 
	// setup parameters for matrix functions (in terms of optical depth)
    
        double eps=sigma_t*delta_x/abs(mu[j]);
	char_funcs<double> p(eps, sigma_t);        // create matrix functions
       
        Type f0p, f1p, f0m, f1m, em, ep;
        double kdx=lambda;
    
        em = Type (cos(-kdx), sin(-kdx));
        ep = Type (cos(kdx), sin(kdx));
    
        f0p = p.p0()*em/(1.0 - em*exp(-eps));
        f1p = p.p1()/(1.0 - em*exp(-eps));
    
        f0m = p.p1()/(1.0 - ep*exp(-eps));
        f1m = p.p0()*ep/(1.0 - ep*exp(-eps));
    
        //Create T positive matrix.
        T_pos[0] = p.p02()*f0p + p.p00();
        T_pos[2] = p.p01()*ep + p.p02()*f1p;
        T_pos[1] = p.p12()*f0p + p.p10()*em;
        T_pos[3] = p.p12()*f1p + p.p11();
    
        //Create T negative matrix.
        T_neg[0] = p.p12()*f0m*ep + p.p11();
        T_neg[2] = p.p12()*f1m*ep + p.p10()*ep;
        T_neg[1] = p.p02()*f0m + p.p01()*em;
        T_neg[3] = p.p02()*f1m + p.p00();
       
        //Define T_pos, T_neg, P matrices as lapack_matrix type.
        lapack_matrix<Type,external>::type T_pos_mat(T_pos,mat_order,mat_order);
        lapack_matrix<Type,external>::type T_neg_mat(T_neg,mat_order,mat_order);
       
        //Summation of T_positive and T_negative.
        scale(T_pos_mat,weight[j]);
        scale(T_neg_mat,weight[j]);
        add(T_pos_mat,T_positive);
        add(T_neg_mat,T_negative);


    }

    //Summation to create T_total matrix.
    add(T_positive,T_total);
    add(T_negative,T_total);
    scale(T_total,(sigma_s/2.0));

    //Return T_total.  The A matrix in the eigensystem A*x=lam*x.
    return T_total;
}
//---------------------------------------------------------------------------//
//                              end of matrix.cc
//---------------------------------------------------------------------------//
