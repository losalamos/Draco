//----------------------------------*-C++-*----------------------------------//
// matrix.t.hh
// John Gulick
// Wed Aug 25 14:05:45 1999
//
//---------------------------------------------------------------------------//
// @> $Id$
//---------------------------------------------------------------------------//
template <class T>
class char_funcs {

  public:
    char_funcs();
    char_funcs(T e, T x);
    ~char_funcs();
    T p00();
    T p01();
    T p02();
    T p10();
    T p11();
    T p12();
    T p0();
    T p1();

  private:
    T xt;
    T eps, eps2, eps3;

};

template <class T>
char_funcs<T>::char_funcs() { //INVALID CONSTRUCTOR
    cerr << "must instantiate char_func class with sigma_t and epsilon" << endl;
    exit(0);
}

template <class T>
char_funcs<T>::char_funcs(T e, T x) {  //CONSTRUCTOR
    eps = e;
    eps2 = e*eps;
    eps3 = e*eps2;
    xt = x;
   
}

template <class T>
char_funcs<T>::~char_funcs() {  //DESTRUCTOR (set private data to zero)
    xt = 0.0;
    eps = 0.0;
    eps2 = 0.0;
    eps3 = 0.0;
}

template <class T>
T char_funcs<T>::p00() {
    T tmp;
    tmp = -exp(-eps)*(6.0/eps3 + 8.0/eps2 + 2.0/eps) 
	+ 6.0/eps3 + 2.0/eps2 - 3.0/eps + 1.0;
    tmp = tmp/xt;
    return(tmp);   
}

template <class T>
T char_funcs<T>::p01() {
    T tmp;
    tmp = exp(-eps)*(6.0/eps3 + 2.0/eps2) - 6.0/eps3 + 4.0/eps2 - 1.0/eps; 
    tmp = tmp/xt;
    return(tmp);   
}

template <class T>
T char_funcs<T>::p02() {
    T tmp;
    tmp = exp(-eps)*(6.0/eps2 + 2.0/eps) - 6.0/eps2 + 4.0/eps; 
    return(tmp);   
}

template <class T>
T char_funcs<T>::p12() {
    T tmp;
    tmp = -exp(-eps)*(6.0/eps2 + 4.0/eps) + 6.0/eps2 - 2.0/eps; 
    return(tmp);   
}

template <class T>
T char_funcs<T>::p10() {
    T tmp;
    tmp = exp(-eps)*(6.0/eps3 + 10.0/eps2 + 4.0/eps) 
	- 6.0/eps3 - 4.0/eps2 + 3.0/eps;
    tmp = tmp/xt;
    return(tmp);   
}

template <class T>
T char_funcs<T>::p11() {
    T tmp;
    tmp = -exp(-eps)*(6.0/eps3 + 4.0/eps2) + 6.0/eps3 - 2.0/eps2 - 1.0/eps + 1.0;
    tmp = tmp/xt;
    return(tmp);   
}

template <class T>
T char_funcs<T>::p0() {
    return ( (1.0 - exp(-eps))/eps - exp(-eps) );
}

template <class T>
T char_funcs<T>::p1() {
   return ( 1.0 + (exp(-eps) - 1.0)/eps );

}
//---------------------------------------------------------------------------//
//                              end of matrix.t.hh
//---------------------------------------------------------------------------//
