#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#define GPU_ERR_CHK(call)                                                   \
    do {                                                                    \
        cudaError_t stat = (call);                                          \
        if (stat != cudaSuccess) {                                          \
            fprintf(stderr, "CUDA error at %s:%d — %s\n",                   \
                    __FILE__, __LINE__, cudaGetErrorString(stat));          \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define CBLAS_ERR_CHK(call)                                                 \
    do {                                                                    \
        cublasStatus_t status = (call);                                     \
        if (status != CUBLAS_STATUS_SUCCESS) {                              \
            fprintf(stderr, "cuBLAS error at %s:%d — code %d\n",            \
                    __FILE__, __LINE__, (int)status);                       \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)

#define TILE_DIM 16

int check_close_floats(const float* arr1, const float* arr2, int len, float tol) {
    for (int k = 0; k < len; k++)
        if (fabsf(arr1[k] - arr2[k]) > tol) return 0;
    return 1;
}

float get_elapsed_ms(cudaEvent_t start_ev, cudaEvent_t stop_ev) {
    float elapsed = 0.0f;
    cudaEventElapsedTime(&elapsed, start_ev, stop_ev);
    return elapsed;
}

__global__ void kernelNaiveGemm(const float* matA, const float* matB, float* matC,
                                int dimM, int dimN, int dimK) {
    int r = blockIdx.y * blockDim.y + threadIdx.y;
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= dimM || c >= dimN) return;
    
    float acc = 0.0f;
    for (int idx = 0; idx < dimK; idx++)
        acc += matA[r * dimK + idx] * matB[idx * dimN + c];
    matC[r * dimN + c] = acc;
}

void execute_naive_gemm(float* d_matA, float* d_matB, float* d_matC,
                        int dimM, int dimN, int dimK, float* time_out) {
    dim3 blk(TILE_DIM, TILE_DIM);
    dim3 grd((dimN + TILE_DIM - 1) / TILE_DIM, (dimM + TILE_DIM - 1) / TILE_DIM);

    cudaEvent_t ev0, ev1;
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);
    kernelNaiveGemm<<<grd, blk>>>(d_matA, d_matB, d_matC, dimM, dimN, dimK);
    cudaEventRecord(ev1);
    cudaEventSynchronize(ev1);
    
    if (time_out) *time_out = get_elapsed_ms(ev0, ev1);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
}

__global__ void kernelTiledGemm(const float* matA, const float* matB, float* matC,
                                int dimM, int dimN, int dimK) {
    __shared__ float subA[TILE_DIM][TILE_DIM];
    __shared__ float subB[TILE_DIM][TILE_DIM];

    int g_row = blockIdx.y * TILE_DIM + threadIdx.y;
    int g_col = blockIdx.x * TILE_DIM + threadIdx.x;
    float acc = 0.0f;

    for (int step = 0; step < (dimK + TILE_DIM - 1) / TILE_DIM; step++) {
        subA[threadIdx.y][threadIdx.x] = (g_row < dimM && step * TILE_DIM + threadIdx.x < dimK)
                                       ? matA[g_row * dimK + step * TILE_DIM + threadIdx.x] : 0.0f;

        subB[threadIdx.y][threadIdx.x] = (g_col < dimN && step * TILE_DIM + threadIdx.y < dimK)
                                       ? matB[(step * TILE_DIM + threadIdx.y) * dimN + g_col] : 0.0f;

        __syncthreads();

        for (int k = 0; k < TILE_DIM; k++) {
            acc += subA[threadIdx.y][k] * subB[k][threadIdx.x];
        }

        __syncthreads();
    }

    if (g_row < dimM && g_col < dimN)
        matC[g_row * dimN + g_col] = acc;
}

void test_tiled_gemm(int dimM, int dimN, int dimK) {
    size_t sz_A = (size_t)dimM * dimK * sizeof(float);
    size_t sz_B = (size_t)dimK * dimN * sizeof(float);
    size_t sz_C = (size_t)dimM * dimN * sizeof(float);

    float *h_matA = (float*)malloc(sz_A);
    float *h_matB = (float*)malloc(sz_B);
    float *h_matC = (float*)malloc(sz_C);
    float *h_chk  = (float*)malloc(sz_C);

    for (int k = 0; k < dimM * dimK; k++) h_matA[k] = ((float)rand() / RAND_MAX - 0.5f);
    for (int k = 0; k < dimK * dimN; k++) h_matB[k] = ((float)rand() / RAND_MAX - 0.5f);

    for (int r = 0; r < dimM; r++) {
        for (int c = 0; c < dimN; c++) {
            float val = 0.0f;
            for (int k = 0; k < dimK; k++) val += h_matA[r * dimK + k] * h_matB[k * dimN + c];
            h_chk[r * dimN + c] = val;
        }
    }

    float *d_matA, *d_matB, *d_matC;
    GPU_ERR_CHK(cudaMalloc(&d_matA, sz_A));
    GPU_ERR_CHK(cudaMalloc(&d_matB, sz_B));
    GPU_ERR_CHK(cudaMalloc(&d_matC, sz_C));
    GPU_ERR_CHK(cudaMemcpy(d_matA, h_matA, sz_A, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_matB, h_matB, sz_B, cudaMemcpyHostToDevice));

    dim3 blk(TILE_DIM, TILE_DIM);
    dim3 grd((dimN + TILE_DIM - 1) / TILE_DIM, (dimM + TILE_DIM - 1) / TILE_DIM);

    cudaEvent_t ev_start, ev_stop;
    cudaEventCreate(&ev_start); cudaEventCreate(&ev_stop);
    cudaEventRecord(ev_start);
    kernelTiledGemm<<<grd, blk>>>(d_matA, d_matB, d_matC, dimM, dimN, dimK);
    cudaEventRecord(ev_stop);
    cudaEventSynchronize(ev_stop);
    float elapsed_t = get_elapsed_ms(ev_start, ev_stop);

    GPU_ERR_CHK(cudaMemcpy(h_matC, d_matC, sz_C, cudaMemcpyDeviceToHost));
    int valid = check_close_floats(h_matC, h_chk, dimM * dimN, 5e-2f);

    double calc_gflops = 2.0 * dimM * dimN * dimK / (elapsed_t / 1000.0) / 1e9;
    printf("  [TiledMatMul] %dx%d@%dx%d  %.2f ms  %.1f GFLOPS  %s\n",
           dimM, dimK, dimK, dimN, elapsed_t, calc_gflops, valid ? "[PASS]" : "[FAIL]");

    cudaEventDestroy(ev_start); cudaEventDestroy(ev_stop);
    cudaFree(d_matA); cudaFree(d_matB); cudaFree(d_matC);
    free(h_matA); free(h_matB); free(h_matC); free(h_chk);
}

void benchmark_gemm_versions(cublasHandle_t cb_handle) {
    int test_dims[] = {128, 256, 512, 1024};
    int num_tests = sizeof(test_dims) / sizeof(test_dims[0]);

    printf("\n  [GemmBenchmark]\n");
    printf("  %6s  %12s  %12s  %12s  %10s\n",
           "Size", "Naive(ms)", "Tiled(ms)", "cuBLAS(ms)", "cuBLAS GFLOPS");
    printf("  %s\n", "--------------------------------------------------------------");

    for (int t = 0; t < num_tests; t++) {
        int dM = test_dims[t], dN = test_dims[t], dK = test_dims[t];
        size_t mem_sz = (size_t)dM * dN * sizeof(float);

        float *d_matA, *d_matB, *d_matC;
        GPU_ERR_CHK(cudaMalloc(&d_matA, mem_sz));
        GPU_ERR_CHK(cudaMalloc(&d_matB, mem_sz));
        GPU_ERR_CHK(cudaMalloc(&d_matC, mem_sz));
        GPU_ERR_CHK(cudaMemset(d_matA, 0, mem_sz));
        GPU_ERR_CHK(cudaMemset(d_matB, 0, mem_sz));

        float ms_naive = 0.0f, ms_tiled = 0.0f, ms_cublas = 0.0f;
        cudaEvent_t ev_st, ev_en;
        GPU_ERR_CHK(cudaEventCreate(&ev_st));
        GPU_ERR_CHK(cudaEventCreate(&ev_en));

        execute_naive_gemm(d_matA, d_matB, d_matC, dM, dN, dK, &ms_naive);

        dim3 blk(TILE_DIM, TILE_DIM);
        dim3 grd((dN + TILE_DIM - 1) / TILE_DIM, (dM + TILE_DIM - 1) / TILE_DIM);
        GPU_ERR_CHK(cudaEventRecord(ev_st));
        kernelTiledGemm<<<grd, blk>>>(d_matA, d_matB, d_matC, dM, dN, dK);
        GPU_ERR_CHK(cudaEventRecord(ev_en));
        GPU_ERR_CHK(cudaEventSynchronize(ev_en));
        ms_tiled = get_elapsed_ms(ev_st, ev_en);

        float s_alpha = 1.0f, s_beta = 0.0f;
        GPU_ERR_CHK(cudaEventRecord(ev_st));
        CBLAS_ERR_CHK(cublasSgemm(cb_handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                 dN, dM, dK, &s_alpha, d_matB, dN, d_matA, dK, &s_beta, d_matC, dN));
        GPU_ERR_CHK(cudaEventRecord(ev_en));
        GPU_ERR_CHK(cudaEventSynchronize(ev_en));
        ms_cublas = get_elapsed_ms(ev_st, ev_en);

        double calc_gf = 2.0 * dM * dN * dK / (ms_cublas / 1000.0 + 1e-9) / 1e9;
        printf("  %6d  %12.2f  %12.2f  %12.2f  %10.1f\n",
               dM, ms_naive, ms_tiled, ms_cublas, calc_gf);

        cudaEventDestroy(ev_st); cudaEventDestroy(ev_en);
        cudaFree(d_matA); cudaFree(d_matB); cudaFree(d_matC);
    }
}

__global__ void kernelMaxPool2x2(const float* in_tensor, float* out_tensor,
                                 int batch, int chn, int ht, int wd) {
    int out_ht = ht / 2;
    int out_wd = wd / 2;

    int b_idx  = blockIdx.z;
    int c_idx  = blockIdx.y;
    int oh_idx = blockIdx.x * blockDim.y + threadIdx.y;
    int ow_idx = threadIdx.x;

    if (oh_idx >= out_ht || ow_idx >= out_wd || b_idx >= batch || c_idx >= chn) return;

    float max_val = -1e30f;
    for (int off_h = 0; off_h < 2; off_h++) {
        for (int off_w = 0; off_w < 2; off_w++) {
            int ih = oh_idx * 2 + off_h;
            int iw = ow_idx * 2 + off_w;
            int lin_idx = ((b_idx * chn + c_idx) * ht + ih) * wd + iw;
            max_val = fmaxf(max_val, in_tensor[lin_idx]);
        }
    }

    out_tensor[((b_idx * chn + c_idx) * out_ht + oh_idx) * out_wd + ow_idx] = max_val;
}

void hostMaxPool(const float* in_data, float* out_data,
                 int b, int c, int h, int w) {
    int h2 = h / 2, w2 = w / 2;
    for (int n_i = 0; n_i < b; n_i++)
      for (int c_i = 0; c_i < c; c_i++)
        for (int oh = 0; oh < h2; oh++)
          for (int ow = 0; ow < w2; ow++) {
              float m_val = -1e30f;
              for (int dh = 0; dh < 2; dh++)
                for (int dw = 0; dw < 2; dw++) {
                    float v_curr = in_data[((n_i * c + c_i) * h + oh * 2 + dh) * w + ow * 2 + dw];
                    if (v_curr > m_val) m_val = v_curr;
                }
              out_data[((n_i * c + c_i) * h2 + oh) * w2 + ow] = m_val;
          }
}

void test_max_pooling(void) {
    int d_n = 4, d_c = 8, d_h = 16, d_w = 16;
    int d_h2 = d_h / 2, d_w2 = d_w / 2;
    size_t in_sz  = (size_t)d_n * d_c * d_h * d_w * sizeof(float);
    size_t out_sz = (size_t)d_n * d_c * d_h2 * d_w2 * sizeof(float);

    float *h_input  = (float*)malloc(in_sz);
    float *h_output = (float*)malloc(out_sz);
    float *h_check  = (float*)malloc(out_sz);
    for (int k = 0; k < d_n * d_c * d_h * d_w; k++) h_input[k] = (float)rand() / RAND_MAX;
    hostMaxPool(h_input, h_check, d_n, d_c, d_h, d_w);

    float *d_input, *d_output;
    GPU_ERR_CHK(cudaMalloc(&d_input,  in_sz));
    GPU_ERR_CHK(cudaMalloc(&d_output, out_sz));
    GPU_ERR_CHK(cudaMemcpy(d_input, h_input, in_sz, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemset(d_output, 0, out_sz));

    dim3 blk(d_w2, 2);
    dim3 grd((d_h2 + 1) / 2, d_c, d_n);
    kernelMaxPool2x2<<<grd, blk>>>(d_input, d_output, d_n, d_c, d_h, d_w);
    GPU_ERR_CHK(cudaMemcpy(h_output, d_output, out_sz, cudaMemcpyDeviceToHost));

    int valid = check_close_floats(h_output, h_check, d_n * d_c * d_h2 * d_w2, 1e-5f);
    printf("  [MaxPool2x2] (%d,%d,%d,%d)->(%d,%d,%d,%d)  %s\n",
           d_n, d_c, d_h, d_w, d_n, d_c, d_h2, d_w2, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_input); cudaFree(d_output);
    free(h_input); free(h_output); free(h_check);
}

__global__ void kernelBatchNormInfer(const float* in_feat, float* out_feat,
                                     const float* gm, const float* bt,
                                     const float* mu,  const float* sigma_sq,
                                     int b_size, int chns, int spatial, float eps_val) {
    int c_idx  = blockIdx.y;
    int sp_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (sp_idx >= spatial || c_idx >= chns) return;

    for (int b_i = 0; b_i < b_size; b_i++) {
        int l_idx = (b_i * chns + c_idx) * spatial + sp_idx;
        float norm_x = (in_feat[l_idx] - mu[c_idx]) / sqrtf(sigma_sq[c_idx] + eps_val);
        out_feat[l_idx] = gm[c_idx] * norm_x + bt[c_idx];
    }
}

void test_batch_norm(void) {
    int n_d=4, c_d=8, h_d=16, w_d=16;
    int sp_d = h_d * w_d;
    float epsilon = 1e-5f;
    size_t data_sz = (size_t)n_d * c_d * sp_d * sizeof(float);
    size_t param_sz = c_d * sizeof(float);

    float *h_tensor  = (float*)malloc(data_sz);
    float *h_result  = (float*)malloc(data_sz);
    float *h_gamma   = (float*)malloc(param_sz);
    float *h_beta    = (float*)malloc(param_sz);
    float *h_mean    = (float*)malloc(param_sz);
    float *h_var     = (float*)malloc(param_sz);
    float *h_check   = (float*)malloc(data_sz);

    for (int k = 0; k < n_d * c_d * sp_d; k++) h_tensor[k] = ((float)rand() / RAND_MAX - 0.5f) * 2.0f;
    for (int c_i = 0; c_i < c_d; c_i++) {
        h_gamma[c_i] = 1.0f;
        h_beta[c_i]  = 0.0f;
        
        double s1 = 0.0, s2 = 0.0;
        for (int b_i = 0; b_i < n_d; b_i++)
            for (int s_i = 0; s_i < sp_d; s_i++) {
                float tmp_v = h_tensor[(b_i * c_d + c_i) * sp_d + s_i];
                s1 += tmp_v; s2 += tmp_v * tmp_v;
            }
        h_mean[c_i] = (float)(s1 / (n_d * sp_d));
        h_var[c_i]  = (float)(s2 / (n_d * sp_d) - h_mean[c_i] * h_mean[c_i]);
        
        for (int b_i = 0; b_i < n_d; b_i++)
            for (int s_i = 0; s_i < sp_d; s_i++) {
                int f_idx = (b_i * c_d + c_i) * sp_d + s_i;
                float x_n = (h_tensor[f_idx] - h_mean[c_i]) / sqrtf(h_var[c_i] + epsilon);
                h_check[f_idx] = h_gamma[c_i] * x_n + h_beta[c_i];
            }
    }

    float *d_tensor, *d_result, *d_gamma, *d_beta, *d_mean, *d_var;
    GPU_ERR_CHK(cudaMalloc(&d_tensor, data_sz));
    GPU_ERR_CHK(cudaMalloc(&d_result, data_sz));
    GPU_ERR_CHK(cudaMalloc(&d_gamma,  param_sz));
    GPU_ERR_CHK(cudaMalloc(&d_beta,   param_sz));
    GPU_ERR_CHK(cudaMalloc(&d_mean,   param_sz));
    GPU_ERR_CHK(cudaMalloc(&d_var,    param_sz));
    
    GPU_ERR_CHK(cudaMemcpy(d_tensor, h_tensor, data_sz, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_gamma,  h_gamma,  param_sz, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_beta,   h_beta,   param_sz, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_mean,   h_mean,   param_sz, cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_var,    h_var,    param_sz, cudaMemcpyHostToDevice));

    int thr_count = 256;
    dim3 blk(thr_count);
    dim3 grd((sp_d + thr_count - 1) / thr_count, c_d);
    kernelBatchNormInfer<<<grd, blk>>>(d_tensor, d_result, d_gamma, d_beta, d_mean, d_var,
                                       n_d, c_d, sp_d, epsilon);
    GPU_ERR_CHK(cudaMemcpy(h_result, d_result, data_sz, cudaMemcpyDeviceToHost));

    int valid = check_close_floats(h_result, h_check, n_d * c_d * sp_d, 1e-4f);
    printf("  [BatchNorm] (%d,%d,%d,%d)  %s\n", n_d, c_d, h_d, w_d, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_tensor); cudaFree(d_result); cudaFree(d_gamma);
    cudaFree(d_beta); cudaFree(d_mean); cudaFree(d_var);
    free(h_tensor); free(h_result); free(h_gamma); free(h_beta);
    free(h_mean); free(h_var); free(h_check);
}

__global__ void kernelConv2dBasic(const float* in_img, const float* w_filter,
                                  float* out_img,
                                  int bs, int in_c, int in_h, int in_w,
                                  int out_c, int k_h, int k_w,
                                  int p_h, int p_w,
                                  int s_h, int s_w) {
    int o_h_dim = (in_h + 2 * p_h - k_h) / s_h + 1;
    int o_w_dim = (in_w + 2 * p_w - k_w) / s_w + 1;

    int b_id  = blockIdx.z;
    int oc_id = blockIdx.y;
    int oh_id = blockIdx.x * blockDim.y + threadIdx.y;
    int ow_id = threadIdx.x;

    if (oh_id >= o_h_dim || ow_id >= o_w_dim || b_id >= bs || oc_id >= out_c) return;

    float acc = 0.0f;
    for (int ic_id = 0; ic_id < in_c; ic_id++) {
        for (int f_h = 0; f_h < k_h; f_h++) {
            for (int f_w = 0; f_w < k_w; f_w++) {
                int curr_ih = oh_id * s_h - p_h + f_h;
                int curr_iw = ow_id * s_w - p_w + f_w;
                if (curr_ih >= 0 && curr_ih < in_h && curr_iw >= 0 && curr_iw < in_w) {
                    float i_val = in_img[((b_id * in_c + ic_id) * in_h + curr_ih) * in_w + curr_iw];
                    float k_val = w_filter[((oc_id * in_c + ic_id) * k_h + f_h) * k_w + f_w];
                    acc += i_val * k_val;
                }
            }
        }
    }
    out_img[((b_id * out_c + oc_id) * o_h_dim + oh_id) * o_w_dim + ow_id] = acc;
}

void test_direct_conv2d(void) {
    int bs=2, ic=1, ih=8, iw=8, oc=4, kh=3, kw=3;
    int ph=1, pw=1, sh=1, sw=1;
    int oh = (ih + 2 * ph - kh) / sh + 1;
    int ow = (iw + 2 * pw - kw) / sw + 1;

    int num_in  = bs * ic * ih * iw;
    int num_w   = oc * ic * kh * kw;
    int num_out = bs * oc * oh * ow;

    float *h_input  = (float*)calloc(num_in, sizeof(float));
    float *h_weight = (float*)calloc(num_w,  sizeof(float));
    float *h_result = (float*)calloc(num_out, sizeof(float));
    for (int k = 0; k < num_in; k++) h_input[k]  = (float)rand() / RAND_MAX;
    for (int k = 0; k < num_w;  k++) h_weight[k] = (float)rand() / RAND_MAX;

    float *d_input, *d_weight, *d_result;
    GPU_ERR_CHK(cudaMalloc(&d_input,  num_in  * sizeof(float)));
    GPU_ERR_CHK(cudaMalloc(&d_weight, num_w   * sizeof(float)));
    GPU_ERR_CHK(cudaMalloc(&d_result, num_out * sizeof(float)));
    
    GPU_ERR_CHK(cudaMemcpy(d_input,  h_input,  num_in  * sizeof(float), cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemcpy(d_weight, h_weight, num_w   * sizeof(float), cudaMemcpyHostToDevice));
    GPU_ERR_CHK(cudaMemset(d_result, 0, num_out * sizeof(float)));

    dim3 blk(ow, 4);
    dim3 grd((oh + 3) / 4, oc, bs);
    kernelConv2dBasic<<<grd, blk>>>(d_input, d_weight, d_result,
                                    bs, ic, ih, iw,
                                    oc, kh, kw,
                                    ph, pw, sh, sw);
    GPU_ERR_CHK(cudaDeviceSynchronize());
    GPU_ERR_CHK(cudaMemcpy(h_result, d_result, num_out * sizeof(float), cudaMemcpyDeviceToHost));

    float val_sum = 0.0f;
    for (int k = 0; k < num_out; k++) val_sum += h_result[k];
    printf("  [Conv2D] out_h=%d out_w=%d  sum=%.2f  %s\n",
           oh, ow, val_sum, val_sum > 0.0f ? "[PASS]" : "[FAIL]");

    cudaFree(d_input); cudaFree(d_weight); cudaFree(d_result);
    free(h_input); free(h_weight); free(h_result);
}

int main(void) {
    printf("\n========================================================\n");
    printf("  CUDA Exercise: Tiled GEMM & CNN Layers\n");
    printf("========================================================\n");

    cudaDeviceProp dev_p;
    GPU_ERR_CHK(cudaGetDeviceProperties(&dev_p, 0));
    printf("  GPU: %s  Peak TFLOPS (FP32): ~%.0f\n\n",
           dev_p.name,
           2.0 * dev_p.multiProcessorCount * dev_p.maxThreadsPerMultiProcessor * dev_p.clockRate * 1e-9);

    cublasHandle_t cb_hndl;
    CBLAS_ERR_CHK(cublasCreate(&cb_hndl));

    {
        int dm=256, dn=256, dk=256;
        float *d_ma, *d_mb, *d_mc;
        GPU_ERR_CHK(cudaMalloc(&d_ma, dm * dk * sizeof(float)));
        GPU_ERR_CHK(cudaMalloc(&d_mb, dk * dn * sizeof(float)));
        GPU_ERR_CHK(cudaMalloc(&d_mc, dm * dn * sizeof(float)));
        float elapsed_msec;
        execute_naive_gemm(d_ma, d_mb, d_mc, dm, dn, dk, &elapsed_msec);
        double gflops_calc = 2.0 * dm * dn * dk / (elapsed_msec / 1000.0) / 1e9;
        printf("  Naive GEMM %dx%d@%dx%d  %.2f ms  %.1f GFLOPS\n", dm, dk, dk, dn, elapsed_msec, gflops_calc);
        cudaFree(d_ma); cudaFree(d_mb); cudaFree(d_mc);
    }

    test_tiled_gemm(512, 512, 512);
    benchmark_gemm_versions(cb_hndl);
    test_max_pooling();
    test_batch_norm();
    test_direct_conv2d();

    cublasDestroy(cb_hndl);

    printf("\n========================================================\n");
    printf("  Execution Complete!\n");
    printf("========================================================\n\n");
    return 0;
}