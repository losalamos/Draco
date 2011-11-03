extern "C"
{
    __global__ void
    vector_add(double const *A, double const *B, double *C, int const N)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x; 
        // if(i%512==0)
        //     printf("index %d\n",i);
        if (i < N) 
            C[i] = A[i] + B[i];
    }
}
