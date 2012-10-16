struct S
{
    S() {};
    int x, y;
};

union U
{
    U() {};
    int i;
    S s;
};

int main(void)
{
    U u;

    return 0;
}
