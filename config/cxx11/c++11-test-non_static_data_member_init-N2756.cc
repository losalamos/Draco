// Non-static data member initializers
// http://events.visualstudio.com/eng/sessions/details?SessionProfile=1876&TrackProfile=1864
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2756.htm
#include <string>
#include <vector>

class widget
{
    int a = 42;
    std::string b = "xyzzy";
    std::vector<int> c = {1,2,3,4};
  public:
    widget() {}
    explicit widget(int val) : a{val} {}
    widget(int i,int j ): c{i,i,i,i,j,j} {}
};

int main(void)
{
    widget w;
    return 0;
}
