#include <stdio.h>
#include <cuda_runtime.h>

int main() {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);

    printf("Number of CUDA devices: %d\n\n", deviceCount);

    for (int i = 0; i < deviceCount; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);

        printf("Device %d: %s\n", i, prop.name);
        printf("Compute Capability: %d.%d\n", prop.major, prop.minor);
        printf("Total Global Memory: %lu bytes (%.2f GB)\n",
               prop.totalGlobalMem,
               prop.totalGlobalMem / (1024.0 * 1024 * 1024));

        printf("Shared Memory per Block: %lu bytes\n", prop.sharedMemPerBlock);
        printf("Constant Memory: %lu bytes\n", prop.totalConstMem);
        printf("Warp Size: %d\n", prop.warpSize);

        printf("Max Threads per Block: %d\n", prop.maxThreadsPerBlock);
        printf("Max Block Dimensions: (%d, %d, %d)\n",
               prop.maxThreadsDim[0],
               prop.maxThreadsDim[1],
               prop.maxThreadsDim[2]);

        printf("Max Grid Dimensions: (%d, %d, %d)\n",
               prop.maxGridSize[0],
               prop.maxGridSize[1],
               prop.maxGridSize[2]);

        printf("Registers per Block: %d\n", prop.regsPerBlock);

        printf("Multiprocessor Count: %d\n", prop.multiProcessorCount);

        printf("Double Precision Support: %s\n",
               (prop.major >= 1 && prop.minor >= 3) ? "Yes" : "No");

        printf("---------------------------------------------------\n");
    }

    return 0;
}