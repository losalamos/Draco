//----------------------------------*-C++-*----------------------------------//
// fourier.cc
// John Gulick
// Tue Aug  3 13:44:08 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "fourier.t.hh"


//***************************************
//Input member function of class fourier.
//***************************************
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
count = 0;
for(int i=0; i <num_delta_x; i++)
    {
    for(int j=0; j<num_points[i]; j++)
        {
          count += 1;
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
spec_rad_dx  = new double[count];
int info;
lambda = 0.;
count  = 0;
max_spec_rad        = 0.;
max_spec_rad_lambda = 0.;

if(analysis_type == "wave")
   {
   num_delta_x = 1;
   num_points[0]=1;
   }

for(int m=0; m<num_delta_x; m++)
    {
    for(int n=0; n<num_points[m]; n++)
        {
	//Set up data for lambda loop.
        count += 1;
        max_spec_rad = 0.;
	max_spec_rad_lambda = 0.;
        lambda = 0.;

        //Start loop over lambda.
        for(int i=0; i<=num_lambda; i++)
            {
            //Zero out matrices and other loop data.
            mtl::set(T_total,0.);
            mtl::set(eigen_vector_left,0.);
            mtl::set(eigen_vector_right,0.);
            mtl::set(eigen_value,0.);
            info = 0.;

            //Get T_total (A in the system A*x=lam*x) based on type of run.
            T_total = matrix_ptr->SCB_matrix(sigma_t, sigma_s, quad_order, lambda, 
					     delta_x);

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
            if(spec_rad_lam[i] > max_spec_rad)
               {
               max_spec_rad        = spec_rad_lam[i];
               max_spec_rad_lambda = lambda;
               }

            //Print out the spectral radius vs lambda data.
            cout <<i <<" " <<lambda <<" " <<spec_rad_lam[i] <<endl;

	    //Update lambda for next loop.
            lambda += lambda_step_size;
	    }

	//Set the spectral radius vs. dx equal to the maximum from lambda
	//loop. 
        spec_rad_dx[count] = max_spec_rad;

	//Print out the spectral radius vs. lambda data.
        cout <<count <<" " <<delta_x <<" " <<spec_rad_dx[count] <<endl;

	//Update delta_x for the next loop.
        delta_x += dx_step_size[m];
	}
    }
}


//---------------------------------------------------------------------------//
//                              end of fourier.cc
//---------------------------------------------------------------------------//
