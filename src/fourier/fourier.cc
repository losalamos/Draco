//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file fourier.cc
 * \author John Gulick
 * \date Tue Aug  3 13:44:08 1999
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "fourier.hh"

//***************************************
//Input member function of class fourier.
//***************************************

//---------------------------------------------------------------------------//
/* 
 * \fn void Fourier::input()
 *
 * Inputs blah blah blah...
 *
 */
//---------------------------------------------------------------------------//

void Fourier::input()
{
vector<string> input_deck;
ifstream infile("input.dat");
istream_iterator<string> ifile(infile);
istream_iterator<string> eos;
copy(ifile, eos, back_inserter(input_deck));
if(input_deck.empty())
   {
   cout << "Error: No data read from input." <<endl;
   }
method            = input_deck[0].c_str();
acceleration_type = input_deck[1].c_str();
analysis_type     = input_deck[2].c_str();
sigma_t           = atof(input_deck[3].c_str());
sigma_s           = atof(input_deck[4].c_str());
delta_x           = atof(input_deck[5].c_str());
quad_order        = atof(input_deck[6].c_str());
lambda_begin      = atof(input_deck[7].c_str());
lambda_end        = atof(input_deck[8].c_str());
lambda_step_size  = atof(input_deck[9].c_str());
delta_x_begin     = atof(input_deck[10].c_str());
num_delta_x       = atoi(input_deck[11].c_str());
region_number = new int[num_delta_x+1];
num_points    = new int[num_delta_x+1];
dx_step_size  = new double[num_delta_x+1];
int count = 0;
 for(int i=12; i<=(12+(3*num_delta_x)); i+=3)
     {
     region_number[count] = atoi(input_deck[i].c_str());
     num_points[count]    = atoi(input_deck[i+1].c_str());
     dx_step_size[count]  = atof(input_deck[i+2].c_str());
     count += 1;
     }
}


//***************************************
//Solve member function of class fourier.
//***************************************
void Fourier::solve()
{
    //Define the complex type and name it Type.
    typedef complex<double> Type;

    //Instantiate a pointer to a method of type "method".
    int mat_order;
    mat_order = 2;
    matrix_ptr = new Matrix(mat_order);
      
    //Determine the number of waves to loop over.
    num_lambda = (int)((lambda_end-lambda_begin)/lambda_step_size);

    //Determine the number of delta_x points to loop over.
    dx_count = 0;
    for(int i=0; i <num_delta_x; i++)
    {
	for(int j=0; j<num_points[i]; j++)
        {
	    dx_count += 1;
        }
    }

    //Define (zero) basic objects for loop.
    spec_rad = new double[num_lambda+1];
    lapack_matrix<Type>::type T_total(mat_order,mat_order); 
    lapack_matrix<Type>::type eigen_vector_left(mat_order,mat_order);
    lapack_matrix<Type>::type eigen_vector_right(mat_order,mat_order);
    mtl::dense1D< Type > eigen_value(mat_order);
    eigen_test   = new double[mat_order];
    spec_rad_lam = new double[num_lambda];
    spec_rad_dx  = new double[dx_count];
    int info;
    lambda = 0.;
    dx_count  = 0;
    max_lam_spec_rad    = 0.;
    max_dx_spec_rad     = 0.;
    max_spec_rad_lambda = 0.;
    max_spec_rad_dx     = 0.;

    if(analysis_type == "wave")
    {
	num_delta_x = 1;
	num_points[0]=1;
    }
    else
    {
	delta_x = dx_step_size[0];
    }

    for(int m=0; m<num_delta_x; m++)
    {
	for(int n=0; n<num_points[m]; n++)
        {
	    //Set up data for lambda loop.
	    dx_count            += 1;
	    max_lam_spec_rad    = 0.;
	    max_spec_rad_lambda = 0.;
	    lambda              = 0.;

	    //Start loop over lambda.
	    if(acceleration_type == "SCB_M4S_DSA")
	    {
		lambda += lambda_step_size;
	    }

	    for(int i=0; i<=num_lambda; i++)
            {
		//Zero out matrices and other loop data.
		mtl::set(T_total,0.);
		mtl::set(eigen_vector_left,0.);
		mtl::set(eigen_vector_right,0.);
		mtl::set(eigen_value,0.);
		info = 0.;

		//Get T_total (A in the system A*x=lam*x) based on type of run.
		//Needs to be done differently...kinda kludgey!
		if(method == "SCB" || method == "LLD")
		{
		    if(acceleration_type == "SI")
		    {
			T_total = matrix_ptr->SCB_matrix(sigma_t, sigma_s,
							 quad_order, lambda,
							 delta_x);
		    }
		    else if(acceleration_type == "SCB_M4S_DSA")
		    {
			T_total = matrix_ptr->SCB_M4S_matrix(sigma_t, sigma_s,
							     quad_order, lambda, 
							     delta_x);
		    }
		}
		else if(method == "UCB")
		{
		    if(acceleration_type == "SI")
		    {
			T_total = matrix_ptr->UCB_matrix(sigma_t, sigma_s,
							 quad_order, lambda,
							 delta_x);
		    }
		    else if(acceleration_type == "SCB_M4S_DSA")
		    {
			T_total = matrix_ptr->UCB_M4S_matrix(sigma_t, sigma_s,
							     quad_order, lambda, 
							     delta_x);
		    }
		}
		else if(method == "LIN_CHAR")
		{
		    T_total = matrix_ptr->LIN_CHAR_matrix(sigma_t, sigma_s,
							  quad_order, lambda,
							  delta_x);
		}


		//Solve for Eigenvalues/Eigenvectors.
		info = geev(GEEV_CALC_BOTH, T_total, eigen_value,eigen_vector_left, 
			    eigen_vector_right);
		if(info > 0)
		{
		    cerr << "Eigensolve failed at geev." << info <<endl;
		}

		//Determine the spectral radius vs. lambda.
		for(int k=0; k<mat_order; k++)
                {
		    eigen_test[k] = abs(eigen_value[k]);
                }
		spec_rad_lam[i] = eigen_test[0];
		for(int l=0; l<mat_order; l++)
                {
		    if(spec_rad_lam[i] < eigen_test[l])
		    {
			spec_rad_lam[i] = eigen_test[l];
		    }
                }
            
		//Pick out largest spectral radius for the lambda loop.
		if(spec_rad_lam[i] > max_lam_spec_rad)
		{
		    max_lam_spec_rad        = spec_rad_lam[i];
		    max_spec_rad_lambda = lambda;
		}

		//Print out the spectral radius vs lambda data.
		//cout <<i <<" " <<lambda <<" " <<spec_rad_lam[i] <<endl;

		//Update lambda for next loop.
		lambda += lambda_step_size;
	    }

	    //Set the spectral radius vs. dx equal to the maximum from lambda
	    //loop. 
	    spec_rad_dx[dx_count] = max_lam_spec_rad;

            if(analysis_type == "delta")
	    {
		//Pick out the largest spectral radius for the delta run.
		if(max_dx_spec_rad < max_lam_spec_rad)
		{
		    max_dx_spec_rad = max_lam_spec_rad;
		    max_spec_rad_dx = delta_x;
		}

		//Print out the spectral radius vs. lambda data.
		//cout <<dx_count <<" " <<delta_x <<" " <<spec_rad_dx[dx_count] 
		//    <<endl;

		//Update delta_x for the next loop.
		delta_x += dx_step_size[m];
	    }
	}
    }
}


//****
//Plot
//****
void Fourier::plot()
{
float *u;
float *v;
for(int i=0; i<2; i++)
    {
    ostringstream  print_option;
    if(i == 0)
       {
       print_option << "/XSERVE";
       }
    else if(i == 1)
       {
       print_option << "/CPS";
       }
    if(analysis_type == "wave")
       {
       u = new float [num_lambda];
       v = new float [num_lambda];
       int step_counter;
       step_counter = 0;
       double lambda_step;
       lambda_step = 0.;
       for(int i=0; i<num_lambda; i++)
           {
	   u[i] = (float) lambda_step;
	   v[i] = (float) spec_rad_lam[i];
	   lambda_step += lambda_step_size;
           step_counter += 1;
	   }

       if(cpgbeg(0,print_option.str().c_str(),1,1) != 1)
          exit(1);
       double scale = 10*(max_lam_spec_rad+0.05);
       int scale_int;
       scale_int = static_cast< int > (scale);
       scale = (scale_int + 1.)/10.;
       cpgenv(lambda_begin, lambda_end, 0.0, scale, 0, 1);
       ostringstream plot_title;
       plot_title <<"Fourier Analysis of " <<method <<" and " <<acceleration_type;
       cpglab("Frequency \\gl", "Spectral Radius \\gr",
	      plot_title.str().c_str());
       cpgsci(11);
       ostringstream plot_info;
       plot_info <<"Scattering Ratio (c) = " 
		 <<(sigma_s/sigma_t)
		 <<", \\gDx = " <<delta_x
	         <<", Quadrature Order = " <<quad_order;
       cpgmtxt("t",-2.5,0.5,0.5,plot_info.str().c_str());
       cpgsci(5);
       ostringstream plot_max;
       plot_max << "Maximum Spectral Radius= "
		<<max_lam_spec_rad <<" at the "
	        <<max_spec_rad_lambda <<" wave.";
       cpgmtxt("t",-1.5,0.5,0.5,plot_max.str().c_str());
       float * max_point_x;
       float * max_point_y;
       max_point_x = new float [1];
       max_point_y = new float [1];
       max_point_x[0] = (float) max_spec_rad_lambda;
       max_point_y[0] = (float) max_lam_spec_rad;
       cpgpt(1,max_point_x,max_point_y,9);
       cpgiden();
       cpgsci(3);
       cpgline(step_counter, u, v);
       cpgend();
       }

    if(analysis_type == "delta")
       {
       u = new float [dx_count];
       v = new float [dx_count];
       delta_x = delta_x_begin;
       int step_counter;
       step_counter = 0;
       for(int m=0; m<num_delta_x; m++)
           {
           for(int n=0; n<num_points[m]; n++)
               {
               u[step_counter] = (float) log(delta_x);
               v[step_counter] = (float) spec_rad_dx[step_counter+1];
               delta_x += dx_step_size[m]; 
               step_counter += 1;
               }
           }
       if(cpgbeg(0,print_option.str().c_str(),1,1) != 1)
          {
          exit(1);
          }  
       double scale = 10*(max_dx_spec_rad+0.05);
       int scale_int;
       scale_int = static_cast< int > (scale);
       scale = (scale_int + 1.)/10.;
       cpgenv(log(delta_x_begin),
	      log(delta_x),0.0,scale,0,10);
       ostringstream plot_title;
       plot_title <<"Fourier Analysis of " <<method <<" and " <<acceleration_type;
       cpglab("Optical Thickness (MFPs)", "Spectral Radius \\gr",
	      plot_title.str().c_str());
       cpgsci(11);
       ostringstream plot_info;
       plot_info <<"Scattering Ratio (c) = " 
		 <<(sigma_s/sigma_t)
	         <<", Quadrature Order = " <<quad_order;

       cpgmtxt("t",-2.5,0.5,0.5,plot_info.str().c_str());
       cpgsci(5);
       ostringstream plot_max;
       plot_max << "Maximum Spectral Radius= " 
		<<max_dx_spec_rad <<" at "
	        <<max_spec_rad_dx <<" mfps.";
       cpgmtxt("t",-1.5,0.5,0.5,plot_max.str().c_str());
       float * max_point_x;
       float * max_point_y;
       max_point_x = new float [1];
       max_point_y = new float [1];
       max_point_x[0] = (float) log(max_spec_rad_dx);
       max_point_y[0] = (float) max_dx_spec_rad;
       cpgpt(1,max_point_x,max_point_y,9);
       cpgiden();
       cpgsci(3);
       cpgline(step_counter,u,v);
       cpgend();
       }
   }
}


//---------------------------------------------------------------------------//
//                              end of fourier.cc
//---------------------------------------------------------------------------//
