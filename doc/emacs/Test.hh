//---------------------------------------------------------------------------//
// Test of EMACS C++ Indentation Style
//---------------------------------------------------------------------------//

namespace rtt_test {

class Test
{
private:
    // some data
    int x;
    double y;
    
public:
    // constructor
    Test();
    
    // functions
    int get_x() { return x; }
    
    // inline functions
    template<class T> inline T do_something(T &);
};

// inline function
template<class T> T Test::do_something(T &x)
{
    double sum;
    for (int i = 0; i < x; i++)
    {
	x[i] += y;
	sum += x[i];
    }
    x.push_back(sum);
    
    // return the array-type
    return x;
}

} // end of namespace rtt_test
