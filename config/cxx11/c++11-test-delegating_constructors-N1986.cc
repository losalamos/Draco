// Allow one constructor to delegate to another
// We can directly initialize bases and members without code duplication.
#include <string>
#include <vector>
#include <iostream>

//---------------------------------------------------------------------------//
class widget
{
    std::string name;
    std::vector<int> data;
  public:
    // ctor 1
    widget( std::string const &s, int size )
        : name(s),
          data(size)
    { /* empty */ }

    // ctor 2
    widget( char const * psz )
        : widget( psz, 100 )
    { /* empty */ }

    // ctor 3
    widget()
        : widget( "STL is great", 42 )
    {
        std::cout << "default" << std::endl;
    }
};
    

//---------------------------------------------------------------------------//
int main(void)
{
    std::string s1("ctor1");
    char const * s2 = s1.c_str();
    widget w1(s1,1);
    widget w2(s2);
    widget w3;
    return 0;
}
