//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstArray.cc
 * \author Giovanni Bavestrelli
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  Array unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../Assert.hh"
#include "ds_test.hh"
#include "../Release.hh"
#include "../ArraySizes.hh"
#include "../Array.hh"

using namespace std;

// function prototype.
void array_tests();

//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << endl
		 << "Version : " << rtt_dsxx::release() << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	cout << "Testing " << argv[0] << ": version " 
	     << rtt_dsxx::release() << endl; 
	array_tests();
    }
    catch( std::exception &err )
    {
	cout << "ERROR: While testing tstArray, " 
	     << err.what()
	     << endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstArray, " 
		  << "An unknown exception was thrown."
		  << std::endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if( rtt_ds_test::passed ) 
    {
        cout << "**** tstSP Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstSP." << endl;

    return 0;
}

//---------------------------------------------------------------------------//

//! \brief The actual tests for the array class.
void array_tests()
{
    using rtt_dsxx::Array;
    using rtt_dsxx::ArraySizes;

   unsigned int k=0,x,y,z;

   // Array sizes 
   // unsigned int Sizes1[]={10};
   unsigned int Sizes2[]={10,20};
   unsigned int Sizes3[]={10,20,30};
   unsigned int Sizes5[]={5,6,7,8,9};

   // Define some arrays 
   Array<int, 2> A2;               // Two-dimensional
   Array<int, 3> A3(Sizes3);       // Three-dimensional
   Array<int, 4> A4(ArraySizes(10)(20)(40)(50));   // Four-dimensional
   const Array<int, 5> A5(Sizes5); // Five-dimensional constant array

   // Traverse the Array a'la STL and fill it in
   for( Array<int, 3>::iterator it=A3.begin(); it<A3.end(); it++ )
       *it=++k;

   // Bounds checking on access.  Require must be enabled to catch
   // out-of-bounds error. (DBC && 1) must be true.
   
   if( DBC && 1 )
   {
       try
       {
	   cout << A3[3][100][1] << endl;
	       FAILMSG("Failed to catch out-of-bounds access!");
       }
       catch ( rtt_dsxx::assertion const & err )
       {
	   // cout << err.what() << endl;
	   PASSMSG("Caught out of bounds access!");
       }
   }
   else
   {
       PASSMSG("out-of-bounds access not tested when Require macro disabled.");
   }

   // Test the dimensions command:
   {
       size_t ndims = A3.dimensions();
       if( ndims == 3 ) 
       { PASSMSG("member function dimensions() reported the correct value."); }
       else
       { FAILMSG("member function dimensions() reported an incorrect value."); }
   }

   // Create some more arrays

   // Test copy constructor
   Array<int,3> CopyOfA3(A3);
   if( CopyOfA3 == A3 )
   { PASSMSG("Array copy constructor works."); }
   else
   { FAILMSG("Array copy constructor fails."); }

   // Test assignment operator
   CopyOfA3=A3;
   if( CopyOfA3 == A3 )
   { PASSMSG("Array assignment operator works."); }
   else
   { FAILMSG("Array assignment operator fails."); }

   // Assignment to self 
   CopyOfA3=CopyOfA3;
   if( CopyOfA3 == A3 )
   { PASSMSG("Assignment to self works."); }
   else
   { FAILMSG("Assignment to self fails."); }

   // Test Swap
   CopyOfA3.swap(A3);
   if( CopyOfA3 == A3 )
   { PASSMSG("Array1.swap(Array2) works."); }
   else
   { FAILMSG("Array1.swap(Array2) fails."); }

   // Empty array should have zero dimensions
   if( A2.size(1) == 0 && A2.size(2) == 0 )
   { PASSMSG("Array.size(int) of empty array works."); }
   else
   { FAILMSG("Array.size(int) of empty array fails."); }  

   // Resize currently empty Array   
   A2.resize( Sizes2 ); 
   if( A2.size(1) == 10 && A2.size(2) == 20 )
   { PASSMSG("Resize of empty array works."); }
   else
   { FAILMSG("Resize of empty array fails."); }  

   // Resize Array, loose elements 
   A3.resize( ArraySizes(10)(20)(30) ); 
   if( A3 != CopyOfA3 )
   { PASSMSG("Default resize command is destructive!"); }
   else
   { FAILMSG("Default resize did not reset the data!"); }  

//---------------------------------------------------------------------------//
// Test indexing
//---------------------------------------------------------------------------//

   // A2 is 10x20 and all zero.
   // A3 is 10x20x30 and each element's value is equal to it's linear index.
   // A4 should be empty.

   A2[1][2]=A2[2][1]+1;         // Indexing 2D Array
   if( A2[1][2] == A2[2][1]+1 )
   { PASSMSG( "Bracket operator works for 2D array." ); }
   else
   { FAILMSG( "Bracket operator fails for 2D array." ); }

   A3[0][0][0]=10;              // Indexing 3D Array
   if( A3[0][0][0] == 10 )
   { PASSMSG( "Bracket operator works for 3D array." ); }
   else
   { FAILMSG( "Bracket operator fails for 3D array." ); }
   
   int old = A4[1][2][3][4]++;      // Indexing 4D Array
   if( A4[1][2][3][4] == old+1 )
   { PASSMSG( "Bracket operator works for 4D array." ); }
   else
   { FAILMSG( "Bracket operator fails for 4D array." ); }

   int aaa=A5[1][2][3][4][5];   // Indexing 5D Array
   if( A5[1][2][3][4][5] == aaa )
   { PASSMSG( "Bracket operator works for 5D array." ); }
   else
   { FAILMSG( "Bracket operator fails for 5D array." ); } 


   k=0;
   bool elem_access_pass(true);
   // Traverse the Array with nested loops
   for (x=0;x<A3.size(1);x++)
    for (y=0;y<A3.size(2);y++)
     for (z=0;z<A3.size(3);z++)
     {
       A3[x][y][z]=++k;

       // Assert that values are the same as when we used iterators above
       if( A3[x][y][z] != CopyOfA3[x][y][z] )
	   elem_access_pass = false;
     } 

   if( elem_access_pass )
   {   PASSMSG("iterator access for element data is good."); }
   else
   {   FAILMSG("iterator access for element data fails."); }

   // Does resize preserve array data?
   unsigned int Sizes3Big[]={20,30,40};
   CopyOfA3.resize(Sizes3Big,0,true);
   CopyOfA3.resize(Sizes3,0,true);
   if( A3 == CopyOfA3 )
   {   PASSMSG("Array.resize() correctly preserves content."); }
   else
   {   PASSMSG("Array.resize() fails to preserve content."); }

   // Call getsubarray and equality for subarrays
   if( A3[0] == CopyOfA3[0] )
   {  PASSMSG("Equality test for 2D SubArray works."); }
   else
   {  FAILMSG("Equality test for 2D SubArray fails."); }

   if( A3[0][0] != CopyOfA3[0][1] ) 
   {  PASSMSG("Inequality test for 1D SubArray works."); }
   else
   {  FAILMSG("Inequality test for 1D SubArray fails."); }

   if( A3[0][0][0] == CopyOfA3[0][0][0] )
   {  PASSMSG("Equality test for 0D SubArray works."); }
   else
   {  FAILMSG("Equality test for 0D SubArray fails."); }

   // Test equality and inequality operators
   old=A3[1][2][3];
   A3[1][2][3]=56;
   if( A3 != CopyOfA3 )
   {  PASSMSG("Inequality operator for Array works."); }
   else
   {  FAILMSG("Inequality operator for Array fails."); } 

   A3[1][2][3]=old;
   if( A3 == CopyOfA3 )
   {  PASSMSG("Equality operator for Array works."); }
   else
   {  FAILMSG("Equality operator for Array fails."); }

   k=0;
   elem_access_pass = true;
   // Traverse Array with nested loops in a much faster way
   for (x=0;x<A3.size(1);x++)
   {
     rtt_dsxx::RefArray<int,2> Z2=A3[x];
     for (y=0;y<A3.size(2);y++)
     {
       rtt_dsxx::RefArray<int,1> Z1=Z2[y];
       for (z=0;z<A3.size(3);z++)
       {
          Z1[z]=++k;

          // Assert that values are the same as when we used iterators above
          if( Z1[z] != A3[x][y][z] || Z1[z] != CopyOfA3[x][y][z] )
	      elem_access_pass = false;
       }
     }
   }   
   if( elem_access_pass )
   {   PASSMSG("iterator access for element data is good (RefArray)."); }
   else
   {   FAILMSG("iterator access for element data fails (RefArray)."); }


   // Play some games with indexing
   old=A3[1][2][3];
   A3[1][2][3]=1;
   A3[1][++A3[1][2][3]][3]=old;
   if( A3[1][2][3] == old )
   {  PASSMSG("Bracket operator test passes."); }
   else
   {  FAILMSG("Bracket operator test fails."); }

   // Play with standard C Arrays
   typedef int ARR[20][30];
   ARR * MyArr = new ARR[10];
   
   k=0;
   elem_access_pass = true;
   // Traverse a C array
   for (x=0;x<10;x++)
    for (y=0;y<20;y++)
     for (z=0;z<30;z++)
     {
        MyArr[x][y][z]=++k;

        if( MyArr[x][y][z] != A3[x][y][z])
	    elem_access_pass = false;
     }
   if( elem_access_pass )
   {   PASSMSG("iterator access for element data is good (C-array)."); }
   else
   {   FAILMSG("iterator access for element data fails (C-Array)."); }

   // Finished playing with C array
   delete [] MyArr;
   MyArr=NULL;

   // Call some member functions
   int s =A3.size();
   if( s == 6000 )
   { PASSMSG("operator size() works."); }
   else
   { FAILMSG("operator size() fails."); }

   // Use STL non mutating algorithm on entire array
   int * pMaximum=std::max_element(A3.begin(),A3.end());
   if( *pMaximum == 6000 )
   { PASSMSG("max_element(Array) works"); }
   else
   { FAILMSG("max_element(Array) fails."); }

   // Use STL mutating algorithm on entire array
   std::replace( A3.begin(), A3.end(), 10, 100 );

   // Use STL algorithm on constant Array
   size_t numtens( std::count(A3.begin(),A3.end(),10) );
   size_t numhund( std::count(A3.begin(),A3.end(),100) );
   if( numtens == 0 && numhund == 2 )
   { PASSMSG("count(b,e,v) works"); }
   else
   { FAILMSG("cout(b,e,v) fails."); }

   // Traverse RefArray using iterator for faster access
   for (Array<int,3>::iterator az=A3[0].begin(), zz=A3[0].end();az!=zz;az++)
       *az=1;

   // Check the size of a RefArray
   if( A3[0].size() == 600 )
   {  PASSMSG("Array.size() operator works."); }
   else
   {  FAILMSG("Array.size() operator fails."); }

   // Try RefArray's size function
   rtt_dsxx::RefArray<int,2> Z2=A3[0];
   if( Z2.size() == Z2.size(1)*Z2.size(2) )
   {  PASSMSG("Array.size(int) operator works."); }
   else
   {  FAILMSG("Array.size(int) operator fails."); }

   // Test some GetRefArray functions
   {
       rtt_dsxx::RefArray<int,3> Z3 = A3.GetRefArray();
       Array<int, 3> const ConstA3(A3);
       rtt_dsxx::RefArray<int,3> const CZ3 = ConstA3.GetRefArray();

       if( Z3 == CZ3 )
	   { PASSMSG("Comparison between Array and Cosnt RefArray works."); }
       else
	   { FAILMSG("Comparison between Array and Const RefArray failed."); }

       if( CZ3.dimensions() == 3 )
	   { PASSMSG("Successfully queried const RefArray for dimensions()."); }
       else
	   { FAILMSG("Successfully queried const RefArray for dimensions()."); }

       rtt_dsxx::RefArray<int,1> Z1=A3[1][1];
       if( Z1.dimensions() == 1 )
           { PASSMSG("Successfully queried RefArray<int,1> for dimensions()."); }
       else
	   { FAILMSG("Successfully queried RefArray<int,1> for dimensions()."); }

       // Test equality operator for RefArray.  Z3 and Z3b are not equal
       // because they are different sizes.
       CopyOfA3.resize(Sizes3Big,0,false);
       rtt_dsxx::RefArray<int,3> Z3b = CopyOfA3.GetRefArray();
       if( Z3 != Z3b )
           { PASSMSG("Successfully detected that two RefArrays have different sizes."); }
       else
           { FAILMSG("Failed to detect that two RefArrays have different sizes."); }
       
       // Repeat above test for Arrays
       if ( A3 != CopyOfA3 )
           { PASSMSG("Successfully detected that two Arrays have different sizes."); }
       else
           { FAILMSG("Failed to detect that two Arrays have different sizes."); }
   }

   // More Array copy constructor tests.
   
   {
       // create a 1D array with no size.
       Array<int,2> B2;
       B2.clear();
       // Test the empty() member function.
       if( B2.empty() )
           { PASSMSG("Successfully detected that the Array is empty."); }
       else
           { FAILMSG("Failed to detect that the Array was empty."); }

       // copy should fail because B1 has not size information.
       Array<int,2> CB2( B2 );
   }

   // Explicit clear
   A3.clear();
   if( A3.size() == 0 )
   {  PASSMSG("Array.clear() operator works."); }
   else
   {  FAILMSG("Array.clear() operator fails."); }

   std::cout<<"Done!\r\n";

   return;
}

//---------------------------------------------------------------------------//
//                        end of tstArray.cc
//---------------------------------------------------------------------------//

