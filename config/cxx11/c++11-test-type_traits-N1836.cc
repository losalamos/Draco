#include <type_traits>

class C
{
  public:
    operator int() const {return 1;}
};

int main(void)
{
    bool test = 
        std::is_integral<int>::value &&
        std::is_class<C>::value &&
        std::is_convertible<C, int>::value;
    return test ? 0 : 1;
}
