//----------------------------------*-C++-*----------------------------------//
// DynArray.pt
// Dave Nystrom
// 14 January 1997
//---------------------------------------------------------------------------//

#include "DynArray.t.hh"

using namespace dsxx;

template class DynArray<int>;
template std::ostream& operator<<( std::ostream& os, 
				  const DynArray<int>& d );

template class DynArray<long>;
template std::ostream& operator<<( std::ostream& os, 
				  const DynArray<long>& d );

template class DynArray<char>;
template std::ostream& operator<<( std::ostream& os, 
				  const DynArray<char>& d );

template class DynArray<float>;
template std::ostream& operator<<( std::ostream& os, 
				  const DynArray<float>& d );

template class DynArray<double>;
template std::ostream& operator<<( std::ostream& os, 
	  	                  const DynArray<double>& d );

#include <string>

template class DynArray<std::string>;
template std::ostream& operator<<( std::ostream& os, 
			          const DynArray<std::string>& d );

// instantiate class DynArray<String>
// instantiate function String ostream& operator<<( ostream& os, const DynArray<String>& d )

template class DynArray<DynArray<int> >;
template class DynArray<DynArray<char> >;
