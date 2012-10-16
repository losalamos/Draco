class C
{
  public:
    explicit operator bool() const { return true; }
};

int main(void)
{
    C c;

    if (c) return 0;
    else return 1;
}
