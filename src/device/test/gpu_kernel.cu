
extern "C"
{
    __global__ void
    sum(int *dest, int a, int b)
    {
        // Assuming a single thread, 1x1x1 block, 1x1 grid
        *dest = a + b;
    }
}
