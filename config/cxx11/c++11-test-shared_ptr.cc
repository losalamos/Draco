#include <memory>

int nfoos = 0;

class Foo
{
  private:
    int v;

  public:
    Foo(void)           : v(0)   { nfoos++; }
    explicit Foo(int i) : v(i)   { nfoos++; }
    Foo(const Foo &f)   : v(f.v) { nfoos++; }
    virtual ~Foo(void) {  nfoos--; }
    virtual int vf() { return v; }
    int f(void) { return v+1; }
};

int main()
{
    std::shared_ptr<int> foo(new int(5));
    return 0;
}
