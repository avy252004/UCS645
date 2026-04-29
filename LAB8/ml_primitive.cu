#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
#define DEF_SIZE (1 << 18)

int verify_floats(const float* arr1, const float* arr2, int len, float tol) {
    for (int idx = 0; idx < len; idx++)
        if (fabsf(arr1[idx] - arr2[idx]) > tol) return 0;
    return 1;
}

__global__ void applyRelu(const float* src, float* dst, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) dst[idx] = fmaxf(0.0f, src[idx]);
}

__global__ void computeSoftmax(const float* in_logits, float* out_probs, int rows, int cols) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= rows) return;

    const float* curr_row = in_logits + r * cols;
    float* out_row  = out_probs  + r * cols;

    float max_v = -1e30f;
    for (int c = 0; c < cols; c++) max_v = fmaxf(max_v, curr_row[c]);

    float sum_e = 0.0f;
    for (int c = 0; c < cols; c++) sum_e += expf(curr_row[c] - max_v);

    for (int c = 0; c < cols; c++) out_row[c] = expf(curr_row[c] - max_v) / sum_e;
}

void execute_softmax_test(void) {
    int r_count = 4, c_count = 10;
    float *h_log = (float*)malloc(r_count * c_count * sizeof(float));
    float *h_prob  = (float*)malloc(r_count * c_count * sizeof(float));
    for (int k = 0; k < r_count * c_count; k++) h_log[k] = (float)rand() / RAND_MAX;

    float *d_log, *d_prob;
    GPU_ASSERT(cudaMalloc(&d_log, r_count * c_count * sizeof(float)));
    GPU_ASSERT(cudaMalloc(&d_prob,  r_count * c_count * sizeof(float)));
    GPU_ASSERT(cudaMemcpy(d_log, h_log, r_count * c_count * sizeof(float),
                          cudaMemcpyHostToDevice));

    int blk_dim = BLOCK_SIZE, grid_dim = (r_count + blk_dim - 1) / blk_dim;
    computeSoftmax<<<grid_dim, blk_dim>>>(d_log, d_prob, r_count, c_count);
    GPU_ASSERT(cudaMemcpy(h_prob, d_prob, r_count * c_count * sizeof(float),
                          cudaMemcpyDeviceToHost));

    int valid = 1;
    for (int r = 0; r < r_count; r++) {
        float s_val = 0.0f;
        for (int c = 0; c < c_count; c++) s_val += h_prob[r * c_count + c];
        if (fabsf(s_val - 1.0f) > 1e-5f) { valid = 0; break; }
    }
    printf("  [Softmax] Row sums = 1.0: %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_log); cudaFree(d_prob);
    free(h_log); free(h_prob);
}

__global__ void applySigmoid(const float* src, float* dst, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        dst[idx] = 1.0f / (1.0f + expf(-src[idx]));
    }
}

void run_sigmoid_test(int num_elements) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_src = (float*)malloc(mem_sz);
    float *h_dst = (float*)malloc(mem_sz);
    float *h_chk = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) {
        h_src[k]  = ((float)rand() / RAND_MAX - 0.5f) * 10.0f;
        h_chk[k] = 1.0f / (1.0f + expf(-h_src[k]));
    }

    float *d_src, *d_dst;
    GPU_ASSERT(cudaMalloc(&d_src, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_dst, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_src, h_src, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_dst, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    applySigmoid<<<grid_dim, blk_dim>>>(d_src, d_dst, num_elements);
    GPU_ASSERT(cudaMemcpy(h_dst, d_dst, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_dst, h_chk, num_elements, 1e-5f);
    printf("  [Sigmoid] %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_src); cudaFree(d_dst);
    free(h_src); free(h_dst); free(h_chk);
}

__global__ void applyTanh(const float* src, float* dst, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        dst[idx] = tanhf(src[idx]);
    }
}

void run_tanh_test(int num_elements) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_src = (float*)malloc(mem_sz);
    float *h_dst = (float*)malloc(mem_sz);
    float *h_chk = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) {
        h_src[k]  = ((float)rand() / RAND_MAX - 0.5f) * 6.0f;
        h_chk[k] = tanhf(h_src[k]);
    }

    float *d_src, *d_dst;
    GPU_ASSERT(cudaMalloc(&d_src, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_dst, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_src, h_src, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_dst, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    applyTanh<<<grid_dim, blk_dim>>>(d_src, d_dst, num_elements);
    GPU_ASSERT(cudaMemcpy(h_dst, d_dst, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_dst, h_chk, num_elements, 1e-5f);
    printf("  [Tanh] %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_src); cudaFree(d_dst);
    free(h_src); free(h_dst); free(h_chk);
}

__global__ void applyLeakyRelu(const float* src, float* dst, float alpha_val, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        dst[idx] = fmaxf(src[idx], alpha_val * src[idx]);
    }
}

void run_leaky_relu_test(int num_elements, float alpha_val) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_src = (float*)malloc(mem_sz);
    float *h_dst = (float*)malloc(mem_sz);
    float *h_chk = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) {
        h_src[k]  = ((float)rand() / RAND_MAX - 0.5f) * 4.0f;
        h_chk[k] = h_src[k] > 0.0f ? h_src[k] : alpha_val * h_src[k];
    }

    float *d_src, *d_dst;
    GPU_ASSERT(cudaMalloc(&d_src, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_dst, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_src, h_src, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_dst, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    applyLeakyRelu<<<grid_dim, blk_dim>>>(d_src, d_dst, alpha_val, num_elements);
    GPU_ASSERT(cudaMemcpy(h_dst, d_dst, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_dst, h_chk, num_elements, 1e-5f);
    printf("  [LeakyReLU] alpha=%.2f  %s\n", alpha_val, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_src); cudaFree(d_dst);
    free(h_src); free(h_dst); free(h_chk);
}

__global__ void calcReluGrad(const float* grad_out, const float* src_fwd, float* grad_in, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        grad_in[idx] = (src_fwd[idx] > 0.0f) ? grad_out[idx] : 0.0f;
    }
}

void run_relu_grad_test(int num_elements) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_fwd = (float*)malloc(mem_sz);
    float *h_g_out = (float*)malloc(mem_sz);
    float *h_g_in = (float*)malloc(mem_sz);
    float *h_chk = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) {
        h_fwd[k]    = ((float)rand() / RAND_MAX - 0.5f) * 4.0f;
        h_g_out[k] = (float)rand() / RAND_MAX;
        h_chk[k]  = h_fwd[k] > 0.0f ? h_g_out[k] : 0.0f;
    }

    float *d_fwd, *d_g_out, *d_g_in;
    GPU_ASSERT(cudaMalloc(&d_fwd, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_g_out, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_g_in, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_fwd, h_fwd, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemcpy(d_g_out, h_g_out, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_g_in, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    calcReluGrad<<<grid_dim, blk_dim>>>(d_g_out, d_fwd, d_g_in, num_elements);
    GPU_ASSERT(cudaMemcpy(h_g_in, d_g_in, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_g_in, h_chk, num_elements, 1e-5f);
    printf("  [ReLUBackward] %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_fwd); cudaFree(d_g_out); cudaFree(d_g_in);
    free(h_fwd); free(h_g_out); free(h_g_in); free(h_chk);
}

__global__ void computeBCE(const float* preds, const float* targets, float* loss_out, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        float p_val = fmaxf(fminf(preds[idx], 1.0f - 1e-7f), 1e-7f);
        loss_out[idx] = -(targets[idx] * logf(p_val) + (1.0f - targets[idx]) * logf(1.0f - p_val));
    }
}

void run_bce_test(int num_elements) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_p = (float*)malloc(mem_sz);
    float *h_t = (float*)malloc(mem_sz);
    float *h_l = (float*)malloc(mem_sz);
    float *h_chk = (float*)malloc(mem_sz);
    for (int k = 0; k < num_elements; k++) {
        h_p[k]   = (float)rand() / RAND_MAX;
        h_t[k] = (rand() % 2) ? 1.0f : 0.0f;
        float p_val = fmaxf(fminf(h_p[k], 1.0f - 1e-7f), 1e-7f);
        h_chk[k]    = -(h_t[k] * logf(p_val) + (1.0f - h_t[k]) * logf(1.0f - p_val));
    }

    float *d_p, *d_t, *d_l;
    GPU_ASSERT(cudaMalloc(&d_p, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_t, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_l, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_p, h_p, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemcpy(d_t, h_t, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_l, 0, mem_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    computeBCE<<<grid_dim, blk_dim>>>(d_p, d_t, d_l, num_elements);
    GPU_ASSERT(cudaMemcpy(h_l, d_l, mem_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_l, h_chk, num_elements, 1e-4f);
    printf("  [BCE-Loss] %s\n", valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_p); cudaFree(d_t); cudaFree(d_l);
    free(h_p); free(h_t); free(h_l); free(h_chk);
}

__global__ void computeCCE(const float* in_logits, const int* in_labels, float* loss_out, int rows, int cols) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= rows) return;

    const float* curr_row = in_logits + r * cols;

    float max_v = -1e30f;
    for (int c = 0; c < cols; c++) {
        max_v = fmaxf(max_v, curr_row[c]);
    }

    float sum_e = 0.0f;
    for (int c = 0; c < cols; c++) {
        sum_e += expf(curr_row[c] - max_v);
    }

    int lbl = in_labels[r];
    loss_out[r] = -(curr_row[lbl] - max_v) + logf(sum_e);
}

void run_cce_test(int num_elements, int cls_count) {
    size_t log_sz = (size_t)num_elements * cls_count * sizeof(float);
    size_t lbl_sz = num_elements * sizeof(int);
    size_t loss_sz = num_elements * sizeof(float);

    float *h_log = (float*)malloc(log_sz);
    int *h_lbl = (int*)malloc(lbl_sz);
    float *h_l = (float*)malloc(loss_sz);
    float *h_chk = (float*)malloc(loss_sz);

    for (int n = 0; n < num_elements; n++) {
        float max_val = -1e30f, sum_val = 0.0f;
        for (int c = 0; c < cls_count; c++) {
            h_log[n*cls_count+c] = ((float)rand() / RAND_MAX - 0.5f) * 4.0f;
            if (h_log[n*cls_count+c] > max_val) max_val = h_log[n*cls_count+c];
        }
        for (int c = 0; c < cls_count; c++) sum_val += expf(h_log[n*cls_count+c] - max_val);
        h_lbl[n] = rand() % cls_count;
        h_chk[n] = -(h_log[n*cls_count + h_lbl[n]] - max_val) + logf(sum_val);
    }

    float *d_log; int *d_lbl; float *d_l;
    GPU_ASSERT(cudaMalloc(&d_log, log_sz));
    GPU_ASSERT(cudaMalloc(&d_lbl, lbl_sz));
    GPU_ASSERT(cudaMalloc(&d_l, loss_sz));
    GPU_ASSERT(cudaMemcpy(d_log, h_log, log_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemcpy(d_lbl, h_lbl, lbl_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_l, 0, loss_sz));

    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    computeCCE<<<grid_dim, blk_dim>>>(d_log, d_lbl, d_l, num_elements, cls_count);
    GPU_ASSERT(cudaMemcpy(h_l, d_l, loss_sz, cudaMemcpyDeviceToHost));

    int valid = verify_floats(h_l, h_chk, num_elements, 1e-4f);
    printf("  [CrossEntropy] N=%d C=%d  %s\n", num_elements, cls_count, valid ? "[PASS]" : "[FAIL]");

    cudaFree(d_log); cudaFree(d_lbl); cudaFree(d_l);
    free(h_log); free(h_lbl); free(h_l); free(h_chk);
}

__global__ void fusedAdamStep(float* weights, const float* grads, float* mom, float* vel,
                              float l_rate, float b1, float b2, float epsilon,
                              float b1_pow, float b2_pow, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) {
        mom[idx] = b1 * mom[idx] + (1.0f - b1) * grads[idx];
        vel[idx] = b2 * vel[idx] + (1.0f - b2) * grads[idx] * grads[idx];

        float m_hat = mom[idx] / (1.0f - b1_pow);
        float v_hat = vel[idx] / (1.0f - b2_pow);

        weights[idx] -= l_rate * m_hat / (sqrtf(v_hat) + epsilon);
    }
}

void run_adam_test(int num_elements) {
    size_t mem_sz = num_elements * sizeof(float);
    float *h_w = (float*)malloc(mem_sz);
    float *h_g = (float*)malloc(mem_sz);
    float *h_m = (float*)calloc(num_elements, sizeof(float));
    float *h_v = (float*)calloc(num_elements, sizeof(float));
    for (int k = 0; k < num_elements; k++) {
        h_w[k] = (float)rand() / RAND_MAX;
        h_g[k] = ((float)rand()/RAND_MAX - 0.5f) * 0.01f;
    }

    float *d_w, *d_g, *d_m, *d_v;
    GPU_ASSERT(cudaMalloc(&d_w, mem_sz)); GPU_ASSERT(cudaMalloc(&d_g, mem_sz));
    GPU_ASSERT(cudaMalloc(&d_m, mem_sz)); GPU_ASSERT(cudaMalloc(&d_v, mem_sz));
    GPU_ASSERT(cudaMemcpy(d_w, h_w, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemcpy(d_g, h_g, mem_sz, cudaMemcpyHostToDevice));
    GPU_ASSERT(cudaMemset(d_m, 0, mem_sz));
    GPU_ASSERT(cudaMemset(d_v, 0, mem_sz));

    float lr_val=1e-3f, beta1=0.9f, beta2=0.999f, eps_val=1e-8f;
    int blk_dim = BLOCK_SIZE, grid_dim = (num_elements + blk_dim - 1) / blk_dim;
    int valid = 1;

    for (int step = 1; step <= 5; step++) {
        float b1_p = powf(beta1, step), b2_p = powf(beta2, step);
        fusedAdamStep<<<grid_dim, blk_dim>>>(d_w, d_g, d_m, d_v, lr_val, beta1, beta2, eps_val, b1_p, b2_p, num_elements);

        for (int k = 0; k < num_elements; k++) {
            h_m[k] = beta1 * h_m[k] + (1.0f - beta1) * h_g[k];
            h_v[k] = beta2 * h_v[k] + (1.0f - beta2) * h_g[k] * h_g[k];
            float m_h = h_m[k] / (1.0f - b1_p);
            float v_h = h_v[k] / (1.0f - b2_p);
            h_w[k] -= lr_val * m_h / (sqrtf(v_h) + eps_val);
        }
    }

    float *h_w_gpu = (float*)malloc(mem_sz);
    GPU_ASSERT(cudaMemcpy(h_w_gpu, d_w, mem_sz, cudaMemcpyDeviceToHost));
    if (!verify_floats(h_w_gpu, h_w, num_elements, 1e-5f)) valid = 0;

    printf("  [Adam] 5 steps  %s\n", valid ? "[PASS]" : "[FAIL]");
    cudaFree(d_w); cudaFree(d_g); cudaFree(d_m); cudaFree(d_v);
    free(h_w); free(h_g); free(h_m); free(h_v); free(h_w_gpu);
}

int main(void) {
    printf("\n========================================================\n");
    printf("  CUDA Exercise: ML Primitives\n");
    printf("========================================================\n");

    cudaDeviceProp dev_prop;
    GPU_ASSERT(cudaGetDeviceProperties(&dev_prop, 0));
    printf("  GPU: %s\n\n", dev_prop.name);

    execute_softmax_test();
    run_sigmoid_test(DEF_SIZE);
    run_tanh_test(DEF_SIZE);
    run_leaky_relu_test(DEF_SIZE, 0.01f);
    run_relu_grad_test(DEF_SIZE);
    run_bce_test(DEF_SIZE);
    run_cce_test(512, 10);
    run_adam_test(1 << 16);

    printf("\n========================================================\n");
    printf("  Execution Complete!\n");
    printf("========================================================\n\n");
    return 0;
}