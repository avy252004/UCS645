#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

#define GPU_ASSERT(call)                                                    \
    do {                                                                    \
        cudaError_t stat = (call);                                          \
        if (stat != cudaSuccess) {                                          \
            fprintf(stderr, "CUDA error at %s:%d — %s\n",                   \
                    __FILE__, __LINE__, cudaGetErrorString(stat));          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define BLOCK_SIZE 256
#define DEF_SIZE (1 << 20)

int verify_floats(const float* arr1, const float* arr2, int len, float tol) {
    for (int idx = 0; idx < len; idx++)
        if (fabsf(arr1[idx] - arr2[idx]) > tol) return 0;
    return 1;
}

__global__ void sharedMemTest(const float* src, float* dst, int len) {
    __shared__ float s_buf[256];
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    s_buf[threadIdx.x] = (idx < len) ? src[idx] : 0.0f;
    __syncthreads();
    if (idx < len)
        dst[idx] = s_buf[threadIdx.x] * 2.0f;
}

__global__ void reduceSumTree(const float* src, float* dst, int len) {
    __shared__ float s_arr[256];
    int t_id = threadIdx.x;
    int idx  = blockIdx.x * blockDim.x + t_id;

    s_arr[t_id] = (idx < len) ? src[idx] : 0.0f;
    __syncthreads();

    for (int step = blockDim.x / 2; step > 0; step >>= 1) {
        if (t_id < step)
            s_arr[t_id] += s_arr[t_id + step];
        __syncthreads();
    }

    if (t_id == 0) dst[blockIdx.x] = s_arr[0];
}

float host_reduce_sum(const float* d_vec, int len) {
    int blk_dim = BLOCK_SIZE;
    int grid_dim  = (len + blk_dim - 1) / blk_dim;
    float *d_parts;
    GPU_ASSERT(cudaMalloc(&d_parts, grid_dim * sizeof(float)));

    reduceSumTree<<<grid_dim, blk_dim>>>(d_vec, d_parts, len);

    float *h_parts = (float*)malloc(grid_dim * sizeof(float));
    GPU_ASSERT(cudaMemcpy(h_parts, d_parts, grid_dim * sizeof(float),
                          cudaMemcpyDeviceToHost));

    float final_sum = 0.0f;
    for (int b = 0; b < grid_dim; b++) final_sum += h_parts[b];

    cudaFree(d_parts);
    free(h_parts);
    return final_sum;
}

void execute_reference_reduce(void) {
    int num_elements = DEF_SIZE;
    float *h_vec = (float*)malloc(num_elements * sizeof(float));
    double sum_host = 0.0;
    for (int k = 0; k < num_elements; k++) {
        h_vec[k] = (float)rand() / RAND_MAX;
        sum_host  += h_vec[k];
    }
    float *d_vec;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(float)));
    GPU_ASSERT(cudaMemcpy(d_vec, h_vec, num_elements * sizeof(float), cudaMemcpyHostToDevice));

    float sum_dev = host_reduce_sum(d_vec, num_elements);
    printf("  [Ref-TreeReduce] GPU=%.2f  CPU=%.2f  Match: %s\n",
           sum_dev, (float)sum_host,
           fabsf(sum_dev - (float)sum_host) < 100.0f ? "[PASS]" : "[FAIL]");

    cudaFree(d_vec);
    free(h_vec);
}

__global__ void copyViaShared(const float* src, float* dst, int len) {
    __shared__ float s_buf[256];
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    s_buf[threadIdx.x] = (idx < len) ? src[idx] : 0.0f;
    __syncthreads();

    if (idx < len) dst[idx] = s_buf[threadIdx.x];
}

void run_shared_copy(void) {
    int num_elements = 1 << 16;
    size_t mem_sz = num_elements * sizeof(float);
    float *h_src  = (float*)malloc(mem_sz);
    float *h_dst = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) h_src[k] = (float)k;

    float *d_src, *d_dst;
    GPU_ASSERT(cudaMalloc(&d_src,  mem_sz));
    GPU_ASSERT(cudaMalloc(&d_dst, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_src, h_src, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_dst, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    copyViaShared<<<grid_dim, blk_dim>>>(d_src, d_dst, num_elements);
    GPU_ASSERT(cudaMemcpy(h_dst, d_dst, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_dst, h_src, num_elements, 1e-5f);
    printf("  [SharedCopy] %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_src); cudaFree(d_dst);
    free(h_src); free(h_dst);
}

__global__ void findMaxBlock(const float* src, float* dst, int len) {
    __shared__ float s_arr[256];
    int t_id = threadIdx.x;
    int idx  = blockIdx.x * blockDim.x + t_id;

    s_arr[t_id] = (idx < len) ? src[idx] : -1e30f;
    __syncthreads();

    for (int step = blockDim.x / 2; step > 0; step >>= 1) {
        if (t_id < step) {
            s_arr[t_id] = fmaxf(s_arr[t_id], s_arr[t_id + step]);
        }
        __syncthreads();
    }

    if (t_id == 0) dst[blockIdx.x] = s_arr[0];
}

void run_max_reduction(void) {
    int num_elements = 1 << 18;
    float *h_vec = (float*)malloc(num_elements * sizeof(float));
    float max_host = -1e30f;
    for (int k = 0; k < num_elements; k++) {
        h_vec[k] = (float)rand() / RAND_MAX * 100.0f;
        if (h_vec[k] > max_host) max_host = h_vec[k];
    }

    float *d_vec;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(float)));
    GPU_ASSERT(cudaMemcpy(d_vec, h_vec, num_elements * sizeof(float), cudaMemcpyHostToDevice));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    float *d_parts;
    GPU_ASSERT(cudaMalloc(&d_parts, grid_dim * sizeof(float)));
    findMaxBlock<<<grid_dim, blk_dim>>>(d_vec, d_parts, num_elements);

    float *h_parts = (float*)malloc(grid_dim * sizeof(float));
    GPU_ASSERT(cudaMemcpy(h_parts, d_parts, grid_dim * sizeof(float),
                          cudaMemcpyDeviceToHost));
    float max_dev = -1e30f;
    for (int b = 0; b < grid_dim; b++)
        if (h_parts[b] > max_dev) max_dev = h_parts[b];

    int valid = fabsf(max_dev - max_host) < 0.01f;
    printf("  [MaxReduce] GPU=%.4f  CPU=%.4f  %s\n",
           max_dev, max_host, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_vec); cudaFree(d_parts);
    free(h_vec); free(h_parts);
}

__global__ void conflictTestKernel(float* src, float* dst, int jump, int len) {
    __shared__ float s_mem[1024];
    int t_id = threadIdx.x;

    s_mem[t_id * jump % 1024] = (t_id < len) ? src[t_id] : 0.0f;
    __syncthreads();
    if (t_id < len)
        dst[t_id] = s_mem[t_id * jump % 1024] * 2.0f;
}

void run_conflict_test(void) {
    int jumps[] = {1, 2, 4, 8, 16, 32};
    int jump_count = sizeof(jumps) / sizeof(jumps[0]);
    int num_elements = 1024, iters = 5000;

    float *d_vec, *d_res;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(float)));
    GPU_ASSERT(cudaMalloc(&d_res,  num_elements * sizeof(float)));
    GPU_ASSERT(cudaMemset(d_vec, 0, num_elements * sizeof(float)));

    cudaEvent_t ev_start, ev_stop;
    GPU_ASSERT(cudaEventCreate(&ev_start));
    GPU_ASSERT(cudaEventCreate(&ev_stop));

    printf("\n  [BankConflictTest]\n");
    printf("  %8s  %12s\n", "Stride", "Time (us)");
    printf("  %s\n", "------------------------");

    float min_us = -1.0f;
    for (int j = 0; j < jump_count; j++) {
        int jump = jumps[j];

        GPU_ASSERT(cudaEventRecord(ev_start));
        for (int r = 0; r < iters; r++) {
            conflictTestKernel<<<1, num_elements>>>(d_vec, d_res, jump, num_elements);
        }
        GPU_ASSERT(cudaEventRecord(ev_stop));
        GPU_ASSERT(cudaEventSynchronize(ev_stop));

        float ms_elapsed;
        GPU_ASSERT(cudaEventElapsedTime(&ms_elapsed, ev_start, ev_stop));
        float us_avg = ms_elapsed * 1000.0f / iters;

        if (min_us < 0.0f) min_us = us_avg;
        printf("  %8d  %12.2f\n", jump, us_avg);
    }

    cudaFree(d_vec); cudaFree(d_res);
    cudaEventDestroy(ev_start); cudaEventDestroy(ev_stop);
}

__global__ void calcHistogram(const int* src, int* hist_out, int len, int bin_count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        atomicAdd(&hist_out[src[idx]], 1);
    }
}

void run_histogram_test(void) {
    int num_elements = 1 << 18, bin_count = 256;
    int *h_vec  = (int*)malloc(num_elements * sizeof(int));
    int *h_res  = (int*)calloc(bin_count, sizeof(int));
    int *h_chk   = (int*)calloc(bin_count, sizeof(int));

    for (int k = 0; k < num_elements; k++) {
        h_vec[k] = rand() % bin_count;
        h_chk[h_vec[k]]++;
    }

    int *d_vec, *d_res;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(int)));
    GPU_ASSERT(cudaMalloc(&d_res, bin_count * sizeof(int)));
    GPU_ASSERT(cudaMemcpy(d_vec, h_vec, num_elements * sizeof(int), cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_res, 0, bin_count * sizeof(int)));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    calcHistogram<<<grid_dim, blk_dim>>>(d_vec, d_res, num_elements, bin_count);
    GPU_ASSERT(cudaMemcpy(h_res, d_res, bin_count * sizeof(int),
                          cudaMemcpyDeviceToHost));

    int valid = 1;
    for (int b = 0; b < bin_count; b++)
        if (h_res[b] != h_chk[b]) { valid = 0; break; }

    printf("  [Histogram] Size=%d bins=%d  %s\n",
           num_elements, bin_count, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_vec); cudaFree(d_res);
    free(h_vec); free(h_res); free(h_chk);
}

__global__ void reduceWarpLevel(const float* src, float* dst, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    float curr = (idx < len) ? src[idx] : 0.0f;

    for (int span = 16; span > 0; span >>= 1) {
        curr += __shfl_down_sync(0xffffffff, curr, span);
    }

    if (threadIdx.x % 32 == 0) atomicAdd(dst, curr);
}

void run_warp_reduce(void) {
    int num_elements = 32;
    float *h_vec = (float*)malloc(num_elements * sizeof(float));
    float sum_host = 0.0f;
    for (int k = 0; k < num_elements; k++) { h_vec[k] = (float)k; sum_host += h_vec[k]; }

    float *d_vec, *d_res;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(float)));
    GPU_ASSERT(cudaMalloc(&d_res, sizeof(float)));
    GPU_ASSERT(cudaMemcpy(d_vec, h_vec, num_elements * sizeof(float), cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_res, 0, sizeof(float)));

    reduceWarpLevel<<<1, 32>>>(d_vec, d_res, num_elements);

    float sum_dev;
    GPU_ASSERT(cudaMemcpy(&sum_dev, d_res, sizeof(float), cudaMemcpyDeviceToHost));
    int valid = fabsf(sum_dev - sum_host) < 0.01f;
    printf("  [WarpReduce] GPU=%.1f  CPU=%.1f  %s\n",
           sum_dev, sum_host, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_vec); cudaFree(d_res);
    free(h_vec);
}

__global__ void calcHistShared(const int* src, int* hist_out, int len, int bin_count) {
    extern __shared__ int s_hist[];

    for (int b = threadIdx.x; b < bin_count; b += blockDim.x) {
        s_hist[b] = 0;
    }
    __syncthreads();

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        atomicAdd(&s_hist[src[idx]], 1);
    }
    __syncthreads();

    for (int b = threadIdx.x; b < bin_count; b += blockDim.x) {
        atomicAdd(&hist_out[b], s_hist[b]);
    }
}

void run_shared_hist(void) {
    int num_elements = 1 << 20, bin_count = 256;
    int *h_vec  = (int*)malloc(num_elements * sizeof(int));
    int *h_res  = (int*)calloc(bin_count, sizeof(int));
    int *h_chk   = (int*)calloc(bin_count, sizeof(int));
    for (int k = 0; k < num_elements; k++) { h_vec[k] = rand() % bin_count; h_chk[h_vec[k]]++; }

    int *d_vec, *d_res;
    GPU_ASSERT(cudaMalloc(&d_vec, num_elements * sizeof(int)));
    GPU_ASSERT(cudaMalloc(&d_res, bin_count * sizeof(int)));
    GPU_ASSERT(cudaMemcpy(d_vec, h_vec, num_elements * sizeof(int), cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_res, 0, bin_count * sizeof(int)));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    int shared_mem_sz = bin_count * sizeof(int);
    calcHistShared<<<grid_dim, blk_dim, shared_mem_sz>>>(d_vec, d_res, num_elements, bin_count);
    GPU_ASSERT(cudaMemcpy(h_res, d_res, bin_count * sizeof(int),
                          cudaMemcpyDeviceToHost));

    int valid = 1;
    for (int b = 0; b < bin_count; b++)
        if (h_res[b] != h_chk[b]) { valid = 0; break; }

    printf("  [SharedHistogram] Size=%d bins=%d  %s\n",
           num_elements, bin_count, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_vec); cudaFree(d_res);
    free(h_vec); free(h_res); free(h_chk);
}

int main(void) {
    printf("\n========================================================\n");
    printf("  CUDA Exercise: Memory Hierarchy & Shared Mem\n");
    printf("========================================================\n");

    cudaDeviceProp dev_prop;
    GPU_ASSERT(cudaGetDeviceProperties(&dev_prop, 0));
    printf("  GPU: %s  Shared mem/block: %zu KB\n\n",
           dev_prop.name, dev_prop.sharedMemPerBlock / 1024);

    execute_reference_reduce();
    run_shared_copy();
    run_max_reduction();
    run_conflict_test();
    run_histogram_test();
    run_warp_reduce();
    run_shared_hist();

    printf("\n========================================================\n");
    printf("  Execution Complete!\n");
    printf("========================================================\n\n");
    return 0;
}