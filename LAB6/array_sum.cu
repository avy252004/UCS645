#include <iostream>
#include <vector>
#include <cuda_runtime.h>


__global__ void sumKernel(float* d_array, float* d_sum, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        atomicAdd(d_sum, d_array[idx]);
    }
}

int main() {
    int N = 10000;
    size_t size = N * sizeof(float);

   
    std::vector<float> h_array(N, 1.5f); 
    float h_sum = 0.0f;

    float *d_array, *d_sum;
    cudaMalloc((void**)&d_array, size);
    cudaMalloc((void**)&d_sum, sizeof(float));

    cudaMemset(d_sum, 0, sizeof(float));

    cudaMemcpy(d_array, h_array.data(), size, cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    sumKernel<<<blocksPerGrid, threadsPerBlock>>>(d_array, d_sum, N);
    cudaDeviceSynchronize(); 
    cudaMemcpy(&h_sum, d_sum, sizeof(float), cudaMemcpyDeviceToHost);

    std::cout << "Total Sum: " << h_sum << std::endl;
    cudaFree(d_array);
    cudaFree(d_sum);

    return 0;
}