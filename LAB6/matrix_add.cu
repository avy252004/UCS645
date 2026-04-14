#include <iostream>
#include <vector>
#include <cuda_runtime.h>

__global__ void matAddKernel(int* A, int* B, int* C, int width, int height) {
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;

    if (col < width && row < height) {
        int idx = row * width + col;
        C[idx] = A[idx] + B[idx];
    }
}

int main() {
    int width = 1024, height = 1024;
    size_t size = width * height * sizeof(int);

    std::vector<int> h_A(width * height, 1);
    std::vector<int> h_B(width * height, 2);
    std::vector<int> h_C(width * height, 0);

    int *d_A, *d_B, *d_C;
    cudaMalloc(&d_A, size);
    cudaMalloc(&d_B, size);
    cudaMalloc(&d_C, size);

    cudaMemcpy(d_A, h_A.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B.data(), size, cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 blocksPerGrid((width + threadsPerBlock.x - 1) / threadsPerBlock.x,
                       (height + threadsPerBlock.y - 1) / threadsPerBlock.y);

    matAddKernel<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, width, height);

    cudaMemcpy(h_C.data(), d_C, size, cudaMemcpyDeviceToHost);

    std::cout << "Top-left element of result matrix: " << h_C[0] << " (Expected 3)" << std::endl;

    cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
    return 0;
}