//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstArray.cc
 * \author Giovanni Bavestrelli
 * \date   Mon Apr 21 16:00:24 MDT 2003
 * \brief  Array unit test.
 * \note   Copyright (c) 2003-2013 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../ArraySizes.hh"
#include "../Array.hh"

using namespace std;
using namespace rtt_dsxx;

// function prototype.
void array_tests( UnitTest &unitTest );

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
    // function setTestName).
    ScalarUnitTest ut( argc, argv, release );
    try
    {
        array_tests(ut);
        ut.status();
    }
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//

//! \brief The actual tests for the array class.
void array_tests( UnitTest &ut )
{
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
   // out-of-bounds error. (DBC & 1) must be true.
   
   if( (DBC & 1) & !(DBC & 8) )
   {
       try
       {
	   cout << A3[3][100][1] << endl;
	       ut.failure("Failed to catch out-of-bounds access!");
       }
       catch ( assertion const & /* err */ )
       {
	   // cout << err.what() << endl;
	   ut.passes("Caught out of bounds access!");
       }
   }
   else
   {
       ut.passes("out-of-bounds access not tested when Require macro disabled.");
   }

   // Test the dimensions command:
   {
       size_t ndims = A3.dimensions();
       if( ndims == 3 ) 
       { ut.passes("member function dimensions() reported the correct value."); }
       else
       { ut.failure("member function dimensions() reported an incorrect value."); }
   }

   // Create some more arrays

   // Test copy constructor
   Array<int,3> CopyOfA3(A3);
   if( CopyOfA3 == A3 )
   { ut.passes("Array copy constructor works."); }
   else
   { ut.failure("Array copy constructor fails."); }

   // Test assignment operator
   CopyOfA3=A3;
   if( CopyOfA3 == A3 )
   { ut.passes("Array assignment operator works."); }
   else
   { ut.failure("Array assignment operator fails."); }

   // Assignment to self 
   CopyOfA3=CopyOfA3;
   if( CopyOfA3 == A3 )
   { ut.passes("Assignment to self works."); }
   else
   { ut.failure("Assignment to self fails."); }

   // Test Swap
   CopyOfA3.swap(A3);
   if( CopyOfA3 == A3 )
   { ut.passes("Array1.swap(Array2) works."); }
   else
   { ut.failure("Array1.swap(Array2) fails."); }

   // Empty array should have zero dimensions
   if( A2.size(1) == 0 && A2.size(2) == 0 )
   { ut.passes("Array.size(int) of empty array works."); }
   else
   { ut.failure("Array.size(int) of empty array fails."); }  

   // Resize currently empty Array   
   A2.resize( Sizes2 ); 
   if( A2.size(1) == 10 && A2.size(2) == 20 )
   { ut.passes("Resize of empty array works."); }
   else
   { ut.failure("Resize of empty array fails."); }  

   // Resize Array, loose elements 
   A3.resize( ArraySizes(10)(20)(30) ); 
   if( A3 != CopyOfA3 )
   { ut.passes("Default resize command is destructive!"); }
   else
   { ut.failure("Default resize did not reset the data!"); }  

//---------------------------------------------------------------------------//
// Test indexing
//---------------------------------------------------------------------------//

   // A2 is 10x20 and all zero.
   // A3 is 10x20x30 and each element's value is equal to it's linear index.
   // A4 should be empty.

   A2[1][2]=A2[2][1]+1;         // Indexing 2D Array
   if( A2[1][2] == A2[2][1]+1 )
   { ut.passes( "Bracket operator works for 2D array." ); }
   else
   { ut.failure( "Bracket operator fails for 2D array." ); }

   A3[0][0][0]=10;              // Indexing 3D Array
   if( A3[0][0][0] == 10 )
   { ut.passes( "Bracket operator works for 3D array." ); }
   else
   { ut.failure( "Bracket operator fails for 3D array." ); }
   
   int old = A4[1][2][3][4]++;      // Indexing 4D Array
   if( A4[1][2][3][4] == old+1 )
   { ut.passes( "Bracket operator works for 4D array." ); }
   else
   { ut.failure( "Bracket operator fails for 4D array." ); }

   int aaa=A5[1][2][3][4][5];   // Indexing 5D Array
   if( A5[1][2][3][4][5] == aaa )
   { ut.passes( "Bracket operator works for 5D array." ); }
   else
   { ut.failure( "Bracket operator fails for 5D array." ); } 


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
   {   ut.passes("iterator access for element data is good."); }
   else
   {   ut.failure("iterator access for element data fails."); }

   // Does resize preserve array data?
   unsigned int Sizes3Big[]={20,30,40};
   CopyOfA3.resize(Sizes3Big,0,true);
   CopyOfA3.resize(Sizes3,0,true);
   if( A3 == CopyOfA3 )
   {   ut.passes("Array.resize() correctly preserves content."); }
   else
   {   ut.passes("Array.resize() fails to preserve content."); }

   // Call getsubarray and equality for subarrays
   if( A3[0] == CopyOfA3[0] )
   {  ut.passes("Equality test for 2D SubArray works."); }
   else
   {  ut.failure("Equality test for 2D SubArray fails."); }

   if( A3[0][0] != CopyOfA3[0][1] ) 
   {  ut.passes("Inequality test for 1D SubArray works."); }
   else
   {  ut.failure("Inequality test for 1D SubArray fails."); }

   if( A3[0][0][0] == CopyOfA3[0][0][0] )
   {  ut.passes("Equality test for 0D SubArray works."); }
   else
   {  ut.failure("Equality test for 0D SubArray fails."); }

   // Test equality and inequality operators
   old=A3[1][2][3];
   A3[1][2][3]=56;
   if( A3 != CopyOfA3 )
   {  ut.passes("Inequality operator for Array works."); }
   else
   {  ut.failure("Inequality operator for Array fails."); } 

   A3[1][2][3]=old;
   if( A3 == CopyOfA3 )
   {  ut.passes("Equality operator for Array works."); }
   else
   {  ut.failure("Equality operator for Array fails."); }

   k=0;
   elem_access_pass = true;
   // Traverse Array with nested loops in a much faster way
   for (x=0;x<A3.size(1);x++)
   {
     RefArray<int,2> Z2=A3[x];
     for (y=0;y<A3.size(2);y++)
     {
       RefArray<int,1> Z1=Z2[y];
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
   {   ut.passes("iterator access for element data is good (RefArray)."); }
   else
   {   ut.failure("iterator access for element data fails (RefArray)."); }


   // Play some games with indexing
   old=A3[1][2][3];
   A3[1][2][3]=1;
   A3[1][++A3[1][2][3]][3]=old;
   if( A3[1][2][3] == old )
   {  ut.passes("Bracket operator test passes."); }
   else
   {  ut.failure("Bracket operator test fails."); }

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
   {   ut.passes("iterator access for element data is good (C-array)."); }
   else
   {   ut.failure("iterator access for element data fails (C-Array)."); }

   // Finished playing with C array
   delete [] MyArr;
   MyArr=NULL;

   // Call some member functions
   int s =A3.size();
   if( s == 6000 )
   { ut.passes("operator size() works."); }
   else
   { ut.failure("operator size() fails."); }

   // Use STL non mutating algorithm on entire array
   int * pMaximum=std::max_element(A3.begin(),A3.end());
   if( *pMaximum == 6000 )
   { ut.passes("max_element(Array) works"); }
   else
   { ut.failure("max_element(Array) fails."); }

   // Use STL mutating algorithm on entire array
   std::replace( A3.begin(), A3.end(), 10, 100 );

   // Use STL algorithm on constant Array
   size_t numtens( std::count(A3.begin(),A3.end(),10) );
   size_t numhund( std::count(A3.begin(),A3.end(),100) );
   if( numtens == 0 && numhund == 2 )
   { ut.passes("count(b,e,v) works"); }
   else
   { ut.failure("cout(b,e,v) fails."); }

   // Traverse RefArray using iterator for faster access
   for (Array<int,3>::iterator az=A3[0].begin(), zz=A3[0].end();az!=zz;az++)
       *az=1;

   // Check the size of a RefArray
   if( A3[0].size() == 600 )
   {  ut.passes("Array.size() operator works."); }
   else
   {  ut.failure("Array.size() operator fails."); }

   // Try RefArray's size function
   RefArray<int,2> Z2=A3[0];
   if( Z2.size() == Z2.size(1)*Z2.size(2) )
   {  ut.passes("Array.size(int) operator works."); }
   else
   {  ut.failure("Array.size(int) operator fails."); }

   // Test some GetRefArray functions
   {
       RefArray<int,3> Z3 = A3.GetRefArray();
       Array<int, 3> const ConstA3(A3);
       RefArray<int,3> const CZ3 = ConstA3.GetRefArray();

       if( Z3 == CZ3 )
	   { ut.passes("Comparison between Array and Cosnt RefArray works."); }
       else
	   { ut.failure("Comparison between Array and Const RefArray failed."); }

       if( CZ3.dimensions() == 3 )
	   { ut.passes("Successfully queried const RefArray for dimensions()."); }
       else
	   { ut.failure("Successfully queried const RefArray for dimensions()."); }

       RefArray<int,1> Z1=A3[1][1];
       if( Z1.dimensions() == 1 )
           { ut.passes("Successfully queried RefArray<int,1> for dimensions()."); }
       else
	   { ut.failure("Successfully queried RefArray<int,1> for dimensions()."); }

       // Test equality operator for RefArray.  Z3 and Z3b are not equal
       // because they are different sizes.
       CopyOfA3.resize(Sizes3Big,0,false);
       RefArray<int,3> Z3b = CopyOfA3.GetRefArray();
       if( Z3 != Z3b )
           { ut.passes("Successfully detected that two RefArrays have different sizes."); }
       else
           { ut.failure("Failed to detect that two RefArrays have different sizes."); }
       
       // Repeat above test for Arrays
       if ( A3 != CopyOfA3 )
           { ut.passes("Successfully detected that two Arrays have different sizes."); }
       else
           { ut.failure("Failed to detect that two Arrays have different sizes."); }
   }

   // More Array copy constructor tests.
   
   {
       // create a 1D array with no size.
       Array<int,2> B2;
       B2.clear();
       // Test the empty() member function.
       if( B2.empty() )
           { ut.passes("Successfully detected that the Array is empty."); }
       else
           { ut.failure("Failed to detect that the Array was empty."); }

       // copy should fail because B1 has not size information.
       Array<int,2> CB2( B2 );
   }

   // Explicit clear
   A3.clear();
   if( A3.size() == 0 )
   {  ut.passes("Array.clear() operator works."); }
   else
   {  ut.failure("Array.clear() operator fails."); }

   std::cout<<"Done!\r\n";

   return;
}

//---------------------------------------------------------------------------//
//                        end of tstArray.cc
//---------------------------------------------------------------------------//

