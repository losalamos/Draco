extern "C"
{
    __global__ void
    vector_add(float const *A, float const *B, float *C, int const N)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x; 
        // if(i%512==0)
        //     printf("index %d\n",i);
        if (i < N) 
            C[i] = A[i] + B[i];
    }
}
