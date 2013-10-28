
int main()
{
    int array[5] = { 1, 2, 3, 4, 5 };
    for (int& x : array)
        x *= 2;
    return 0;
}
