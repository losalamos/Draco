//----------------------------------*-C++-*----------------------------------//
// ContainerTraits.hh
// Randy M. Roberts
// Wed May  6 14:50:14 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __traits_ContainerTraits_hh__
#define __traits_ContainerTraits_hh__

#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class ContainerTraits - 

// 
//===========================================================================//

template<class T1, class T2=T1>
class ContainerTraits { };

template <>
class ContainerTraits<double>
{
  public:
    typedef const double *const_iterator;
    typedef double *iterator;
    static inline iterator begin(double &a)
    {
	return &a;
    }
    static inline const_iterator begin(const double &a)
    {
	return &a;
    }
    static inline iterator end(double &a)
    {
	return &a+1;
    }
    static inline const_iterator end(const double &a)
    {
	return &a+1;
    }
    static inline bool conformal(double a, double b)
    {
	return true;
    }
};

template <>
class ContainerTraits<int>
{
  public:
    typedef const int *const_iterator;
    typedef int *iterator;
    static inline iterator begin(int &a)
    {
	return &a;
    }
    static inline const_iterator begin(const int &a)
    {
	return &a;
    }
    static inline iterator end(int &a)
    {
	return &a+1;
    }
    static inline const_iterator end(const int &a)
    {
	return &a+1;
    }
    static inline bool conformal(int a, int b)
    {
	return true;
    }
};

template <>
class ContainerTraits<int, double>
{
  public:
    static inline bool conformal(int a, double b)
    {
	return true;
    }
};

template <>
class ContainerTraits<double, int>
{
  public:
    static inline bool conformal(double a, int b)
    {
	return true;
    }
};

template<class T1>
class ContainerTraits<std::vector<T1> >
{
  public:
    typedef typename std::vector<T1>::iterator iterator;
    typedef typename std::vector<T1>::const_iterator const_iterator;
    static inline iterator begin(std::vector<T1> &a)
    {
	return a.begin();
    }
    static inline const_iterator begin(const std::vector<T1> &a)
    {
	return a.begin();
    }
    static inline iterator end(std::vector<T1> &a)
    {
	return a.end();
    }
    static inline const_iterator end(const std::vector<T1> &a)
    {
	return a.end();
    }
    static inline bool conformal(const std::vector<T1> &a,
				 const std::vector<T1> &b)
    {
	return a.size() == b.size();
    }
};

template<class T1, class T2>
class ContainerTraits<std::vector<T1>, std::vector<T2> >
{
  public:
    static inline bool conformal(const std::vector<T1> &a,
				 const std::vector<T2> &b)
    {
	return a.size() == b.size();
    }
};

END_NS_XTM  // namespace XTM

#endif                          // __traits_ContainerTraits_hh__

//---------------------------------------------------------------------------//
//                              end of traits/ContainerTraits.hh
//---------------------------------------------------------------------------//
