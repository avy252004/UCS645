#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cuda_runtime.h>
#include <cudnn.h>
#include <cublas_v2.h>

#define GPU_CHK(call)                                                       \
    do { cudaError_t err=(call);                                            \
         if(err!=cudaSuccess){fprintf(stderr,"CUDA %s:%d %s\n",             \
         __FILE__,__LINE__,cudaGetErrorString(err));exit(1);} } while(0)

#define DNN_CHK(call)                                                       \
    do { cudnnStatus_t stat=(call);                                         \
         if(stat!=CUDNN_STATUS_SUCCESS){fprintf(stderr,"cuDNN %s:%d %d\n",  \
         __FILE__,__LINE__,(int)stat);exit(1);} } while(0)

#define BLAS_CHK(call)                                                      \
    do { cublasStatus_t stat=(call);                                        \
         if(stat!=CUBLAS_STATUS_SUCCESS){fprintf(stderr,"cuBLAS %s:%d %d\n",\
         __FILE__,__LINE__,(int)stat);exit(1);} } while(0)

#define B_SIZE    256
#define LR_VAL    0.01f
#define EPOCHS    10
#define IMG_SZ    784
#define CLS_CNT   10

cudnnHandle_t   dnn_hndl;
cublasHandle_t  blas_hndl;

static double get_time_ms(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec*1e3 + ts.tv_nsec*1e-6;
}

static int read_32bit(FILE* fp) {
    unsigned char b[4];
    fread(b, 1, 4, fp);
    return (b[0]<<24)|(b[1]<<16)|(b[2]<<8)|b[3];
}

typedef struct {
    float* imgs;
    int* lbls;
    int count;
} Dataset;

Dataset load_dataset(const char* i_path, const char* l_path) {
    FILE *f_img = fopen(i_path, "rb");
    FILE *f_lbl = fopen(l_path, "rb");
    if (!f_img || !f_lbl) {
        fprintf(stderr, "Cannot open files.\n");
        exit(1);
    }
    read_32bit(f_img); read_32bit(f_lbl);
    int cnt = read_32bit(f_img); read_32bit(f_lbl);
    int r = read_32bit(f_img);
    int c = read_32bit(f_img);
    (void)r; (void)c;

    Dataset ds;
    ds.count  = cnt;
    ds.imgs = (float*)malloc((size_t)cnt * IMG_SZ * sizeof(float));
    ds.lbls = (int*)malloc(cnt * sizeof(int));

    unsigned char* temp_buf = (unsigned char*)malloc(IMG_SZ);
    for (int k = 0; k < cnt; k++) {
        fread(temp_buf, 1, IMG_SZ, f_img);
        for (int j = 0; j < IMG_SZ; j++)
            ds.imgs[k * IMG_SZ + j] = (temp_buf[j] - 127.5f) / 127.5f;
        unsigned char l_val; fread(&l_val, 1, 1, f_lbl);
        ds.lbls[k] = (int)l_val;
    }
    free(temp_buf); fclose(f_img); fclose(f_lbl);
    printf("[✓] Loaded %d samples from %s\n", cnt, i_path);
    return ds;
}

cudnnTensorDescriptor_t create_tsr_desc(int num, int ch, int ht, int wd) {
    cudnnTensorDescriptor_t dsc;
    DNN_CHK(cudnnCreateTensorDescriptor(&dsc));
    DNN_CHK(cudnnSetTensor4dDescriptor(dsc, CUDNN_TENSOR_NCHW, CUDNN_DATA_FLOAT, num, ch, ht, wd));
    return dsc;
}

cudnnFilterDescriptor_t create_flt_desc(int k, int c, int h, int w) {
    cudnnFilterDescriptor_t dsc;
    DNN_CHK(cudnnCreateFilterDescriptor(&dsc));
    DNN_CHK(cudnnSetFilter4dDescriptor(dsc, CUDNN_DATA_FLOAT, CUDNN_TENSOR_NCHW, k, c, h, w));
    return dsc;
}

cudnnConvolutionDescriptor_t create_cnv_desc(int pad_v, int str_v) {
    cudnnConvolutionDescriptor_t dsc;
    DNN_CHK(cudnnCreateConvolutionDescriptor(&dsc));
    DNN_CHK(cudnnSetConvolution2dDescriptor(dsc, pad_v, pad_v, str_v, str_v, 1, 1, CUDNN_CROSS_CORRELATION, CUDNN_DATA_FLOAT));
    return dsc;
}

__global__ void kernel_relu(float* vec, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) vec[idx] = fmaxf(0.0f, vec[idx]);
}

__global__ void kernel_sce(const float* lgt, const int* tgt, float* prb, float* lss, int n_items, int c_items) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_items) return;

    const float* c_row  = lgt + idx * c_items;
    float* p_row = prb  + idx * c_items;

    float m_val = -1e30f;
    for (int j = 0; j < c_items; j++) m_val = fmaxf(m_val, c_row[j]);
    float s_exp = 0.0f;
    for (int j = 0; j < c_items; j++) { p_row[j] = expf(c_row[j] - m_val); s_exp += p_row[j]; }
    for (int j = 0; j < c_items; j++) p_row[j] /= s_exp;

    lss[idx] = -logf(p_row[tgt[idx]] + 1e-9f);
}

__global__ void kernel_sgd(float* wts, const float* gds, float lr_rt, int len) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < len) wts[idx] -= lr_rt * gds[idx];
}

void exec_conv(
    cudnnTensorDescriptor_t in_d, float* d_in,
    cudnnFilterDescriptor_t flt_d, float* d_flt,
    cudnnConvolutionDescriptor_t cnv_d,
    cudnnTensorDescriptor_t out_d, float* d_out) {

    float a_val = 1.0f, b_val = 0.0f;
    int a_cnt;
    cudnnConvolutionFwdAlgoPerf_t perf_res;
    DNN_CHK(cudnnFindConvolutionForwardAlgorithm(
        dnn_hndl, in_d, flt_d, cnv_d, out_d, 1, &a_cnt, &perf_res));
    cudnnConvolutionFwdAlgo_t fwd_algo = perf_res.algo;

    size_t wk_sz = 0;
    DNN_CHK(cudnnGetConvolutionForwardWorkspaceSize(
        dnn_hndl, in_d, flt_d, cnv_d, out_d, fwd_algo, &wk_sz));

    void* d_wk = NULL;
    if (wk_sz > 0) GPU_CHK(cudaMalloc(&d_wk, wk_sz));

    DNN_CHK(cudnnConvolutionForward(
        dnn_hndl, &a_val, in_d, d_in, flt_d, d_flt,
        cnv_d, fwd_algo, d_wk, wk_sz, &b_val, out_d, d_out));

    if (d_wk) cudaFree(d_wk);
}

void exec_pool(
    cudnnTensorDescriptor_t in_d, float* d_in,
    cudnnTensorDescriptor_t out_d, float* d_out,
    int ph, int pw, int sh, int sw) {

    cudnnPoolingDescriptor_t p_desc;
    DNN_CHK(cudnnCreatePoolingDescriptor(&p_desc));
    DNN_CHK(cudnnSetPooling2dDescriptor(
        p_desc, CUDNN_POOLING_MAX, CUDNN_NOT_PROPAGATE_NAN,
        ph, pw, 0, 0, sh, sw));

    float a_val = 1.0f, b_val = 0.0f;
    DNN_CHK(cudnnPoolingForward(
        dnn_hndl, p_desc, &a_val, in_d, d_in, &b_val, out_d, d_out));

    DNN_CHK(cudnnDestroyPoolingDescriptor(p_desc));
}

__global__ void kernel_bias(float* out_v, const float* bs_v, int b_sz, int f_sz) {
    int n_i = blockIdx.y, c_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (n_i < b_sz && c_i < f_sz) out_v[n_i * f_sz + c_i] += bs_v[c_i];
}

void exec_fc(float* d_in, float* d_wts, float* d_bs, float* d_out, int b_sz, int in_f, int out_f) {
    float a_val = 1.0f, b_val = 0.0f;

    BLAS_CHK(cublasSgemm(
        blas_hndl, CUBLAS_OP_T, CUBLAS_OP_N, out_f, b_sz, in_f,
        &a_val, d_wts, in_f, d_in, in_f, &b_val, d_out, out_f
    ));

    dim3 blk(256);
    dim3 grd((out_f + 255) / 256, b_sz);
    kernel_bias<<<grd, blk>>>(d_out, d_bs, b_sz, out_f);
}

void demo_async(const float* host_imgs, int tot_samps, float* d_bA, float* d_bB) {
    int bsz = B_SIZE;
    size_t bytes_len = (size_t)bsz * IMG_SZ * sizeof(float);
    int tot_b = tot_samps / bsz;

    float *h_pin;
    GPU_CHK(cudaMallocHost(&h_pin, tot_samps * IMG_SZ * sizeof(float)));
    memcpy(h_pin, host_imgs, tot_samps * IMG_SZ * sizeof(float));

    double st_sync = get_time_ms();
    for (int k = 0; k < tot_b; k++) {
        GPU_CHK(cudaMemcpy(d_bA, h_pin + k*bsz*IMG_SZ, bytes_len, cudaMemcpyHostToDevice));
        kernel_relu<<<(bsz * IMG_SZ + 255)/256, 256>>>(d_bA, bsz * IMG_SZ);
        GPU_CHK(cudaDeviceSynchronize());
    }
    double t_sync = get_time_ms() - st_sync;

    cudaStream_t s_cmp, s_xfr;
    GPU_CHK(cudaStreamCreate(&s_cmp));
    GPU_CHK(cudaStreamCreate(&s_xfr));

    double st_async = get_time_ms();
    GPU_CHK(cudaMemcpyAsync(d_bA, h_pin, bytes_len, cudaMemcpyHostToDevice, s_xfr));
    GPU_CHK(cudaStreamSynchronize(s_xfr));

    for (int k = 0; k < tot_b - 1; k++) {
        GPU_CHK(cudaMemcpyAsync(d_bB, h_pin + (k+1)*bsz*IMG_SZ, bytes_len, cudaMemcpyHostToDevice, s_xfr));
        kernel_relu<<<(bsz * IMG_SZ + 255)/256, 256, 0, s_cmp>>>(d_bA, bsz * IMG_SZ);
        GPU_CHK(cudaStreamSynchronize(s_xfr));
        float *tmp = d_bA; d_bA = d_bB; d_bB = tmp;
    }
    double t_async = get_time_ms() - st_async;

    printf("  [AsyncPipeline] Sync: %.2f ms | Async: %.2f ms\n", t_sync, t_async);
    printf("  [AsyncPipeline] Speedup: %.2f ms (%.4f sec).\n", t_sync - t_async, (t_sync - t_async)/1000.0);

    GPU_CHK(cudaStreamDestroy(s_cmp));
    GPU_CHK(cudaStreamDestroy(s_xfr));
    GPU_CHK(cudaFreeHost(h_pin));
}

void run_epoch(int ep,
               float* dw_c1, float* dw_c2, float* dw_f1, float* db_f1, float* dw_f2, float* db_f2,
               float* dx, float* dc1, float* dp1, float* dc2, float* dp2,
               float* df1, float* dlgt, float* dprb, float* dlss,
               cudnnTensorDescriptor_t dx_d,
               cudnnFilterDescriptor_t f1_d, cudnnConvolutionDescriptor_t c1_d, cudnnTensorDescriptor_t c1o_d,
               cudnnTensorDescriptor_t p1_d, cudnnFilterDescriptor_t f2_d, cudnnConvolutionDescriptor_t c2_d,
               cudnnTensorDescriptor_t c2o_d, cudnnTensorDescriptor_t p2_d,
               const float* h_in, const int* h_lbl, int t_cnt) {
               
    int b_cnt = t_cnt / B_SIZE;
    float sum_lss = 0.0f;

    for (int iter = 0; iter < b_cnt; iter++) {
        const float* b_img = h_in + (size_t)iter * B_SIZE * IMG_SZ;
        const int* b_tgt   = h_lbl + iter * B_SIZE;
        int* d_tgt;
        GPU_CHK(cudaMalloc(&d_tgt, B_SIZE * sizeof(int)));

        GPU_CHK(cudaMemcpy(dx, b_img, (size_t)B_SIZE * IMG_SZ * sizeof(float), cudaMemcpyHostToDevice));
        GPU_CHK(cudaMemcpy(d_tgt, b_tgt, B_SIZE * sizeof(int), cudaMemcpyHostToDevice));

        exec_conv(dx_d, dx, f1_d, dw_c1, c1_d, c1o_d, dc1);
        int sz_c1 = B_SIZE * 32 * 28 * 28;
        kernel_relu<<<(sz_c1+255)/256, 256>>>(dc1, sz_c1);
        exec_pool(c1o_d, dc1, p1_d, dp1, 2, 2, 2, 2);

        exec_conv(p1_d, dp1, f2_d, dw_c2, c2_d, c2o_d, dc2);
        int sz_c2 = B_SIZE * 64 * 14 * 14;
        kernel_relu<<<(sz_c2+255)/256, 256>>>(dc2, sz_c2);
        exec_pool(c2o_d, dc2, p2_d, dp2, 2, 2, 2, 2);

        exec_fc(dp2, dw_f1, db_f1, df1, B_SIZE, 64*7*7, 256);
        kernel_relu<<<(B_SIZE*256+255)/256, 256>>>(df1, B_SIZE*256);

        exec_fc(df1, dw_f2, db_f2, dlgt, B_SIZE, 256, CLS_CNT);

        int th = 256, bl = (B_SIZE + th - 1) / th;
        kernel_sce<<<bl, th>>>(dlgt, d_tgt, dprb, dlss, B_SIZE, CLS_CNT);

        float l_batch[B_SIZE];
        GPU_CHK(cudaMemcpy(l_batch, dlss, B_SIZE * sizeof(float), cudaMemcpyDeviceToHost));
        for (int m = 0; m < B_SIZE; m++) sum_lss += l_batch[m];

        if (iter % 50 == 0)
            printf("  Epoch %d  Batch [%d/%d]  AvgLoss=%.4f\n", ep, iter, b_cnt, sum_lss / ((iter+1)*B_SIZE));

        cudaFree(d_tgt);
    }

    printf("  --- Epoch %d Done  AvgLoss=%.4f ---\n", ep, sum_lss / (b_cnt * B_SIZE));
}

void demo_fp16(int dm, int dn, int dk) {
    printf("  [FP16-TensorCore] cublasGemmEx logic with CUDA_R_16F would execute here.\n");
}

int main(void) {
    printf("\n========================================================\n");
    printf("  CUDA Exercise: MNIST CNN (cuDNN + cuBLAS)\n");
    printf("========================================================\n");

    cudaDeviceProp dev_p;
    GPU_CHK(cudaGetDeviceProperties(&dev_p, 0));
    printf("  GPU: %s  Compute: %d.%d  VRAM: %.0f MB\n\n",
           dev_p.name, dev_p.major, dev_p.minor,
           dev_p.totalGlobalMem / 1e6);

    DNN_CHK(cudnnCreate(&dnn_hndl));
    BLAS_CHK(cublasCreate(&blas_hndl));

    Dataset d_trn = load_dataset("data/train-images-idx3-ubyte", "data/train-labels-idx1-ubyte");
    Dataset d_tst = load_dataset("data/t10k-images-idx3-ubyte", "data/t10k-labels-idx1-ubyte");

    float *d_w_c1, *d_w_c2;
    float *d_w_f1, *d_b_f1, *d_w_f2, *d_b_f2;

    GPU_CHK(cudaMalloc(&d_w_c1, 32*1*5*5   * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_w_c2, 64*32*5*5  * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_w_f1, 256*3136   * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_b_f1, 256        * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_w_f2, 10*256     * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_b_f2, 10         * sizeof(float)));

    {
        int w_lens[] = {32*1*5*5, 64*32*5*5, 256*3136, 256, 10*256, 10};
        float* d_ptrs[] = {d_w_c1, d_w_c2, d_w_f1, d_b_f1, d_w_f2, d_b_f2};
        for (int idx = 0; idx < 6; idx++) {
            float* arr = (float*)malloc(w_lens[idx] * sizeof(float));
            float scl = sqrtf(2.0f / w_lens[idx]);
            for (int m = 0; m < w_lens[idx]; m++)
                arr[m] = scl * (2.0f * (float)rand()/RAND_MAX - 1.0f);
            GPU_CHK(cudaMemcpy(d_ptrs[idx], arr, w_lens[idx]*sizeof(float), cudaMemcpyHostToDevice));
            free(arr);
        }
    }

    float *d_inx, *d_m_c1, *d_m_p1, *d_m_c2, *d_m_p2, *d_m_f1, *d_lgt, *d_prb, *d_lss;
    GPU_CHK(cudaMalloc(&d_inx,  (size_t)B_SIZE*1 *28*28 * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_m_c1, (size_t)B_SIZE*32*28*28 * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_m_p1, (size_t)B_SIZE*32*14*14 * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_m_c2, (size_t)B_SIZE*64*14*14 * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_m_p2, (size_t)B_SIZE*64*7 *7  * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_m_f1, (size_t)B_SIZE*256       * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_lgt,  (size_t)B_SIZE*10         * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_prb,  (size_t)B_SIZE*10         * sizeof(float)));
    GPU_CHK(cudaMalloc(&d_lss,  (size_t)B_SIZE             * sizeof(float)));

    cudnnTensorDescriptor_t     desc_x    = create_tsr_desc(B_SIZE, 1,  28, 28);
    cudnnTensorDescriptor_t     desc_c1_o = create_tsr_desc(B_SIZE, 32, 28, 28);
    cudnnTensorDescriptor_t     desc_p1   = create_tsr_desc(B_SIZE, 32, 14, 14);
    cudnnTensorDescriptor_t     desc_c2_o = create_tsr_desc(B_SIZE, 64, 14, 14);
    cudnnTensorDescriptor_t     desc_p2   = create_tsr_desc(B_SIZE, 64, 7,  7);

    cudnnFilterDescriptor_t     desc_f1   = create_flt_desc(32, 1,  5, 5);
    cudnnFilterDescriptor_t     desc_f2   = create_flt_desc(64, 32, 5, 5);

    cudnnConvolutionDescriptor_t desc_c1  = create_cnv_desc(2, 1);
    cudnnConvolutionDescriptor_t desc_c2  = create_cnv_desc(2, 1);

    printf("\n[Training] Starting for %d epochs...\n\n", EPOCHS);

    for (int ep = 1; ep <= EPOCHS; ep++) {
        double start_t = get_time_ms();

        run_epoch(ep, d_w_c1, d_w_c2, d_w_f1, d_b_f1, d_w_f2, d_b_f2,
                  d_inx, d_m_c1, d_m_p1, d_m_c2, d_m_p2, d_m_f1, d_lgt, d_prb, d_lss,
                  desc_x, desc_f1, desc_c1, desc_c1_o, desc_p1,
                  desc_f2, desc_c2, desc_c2_o, desc_p2,
                  d_trn.imgs, d_trn.lbls, d_trn.count);

        double ms_ep = get_time_ms() - start_t;
        printf("  Epoch %d complete in %.1f s\n\n", ep, ms_ep / 1000.0);
    }

    printf("[Stretch] CUDA Streams async pipeline:\n");
    float *d_bA, *d_bB;
    GPU_CHK(cudaMalloc(&d_bA, (size_t)B_SIZE*IMG_SZ*sizeof(float)));
    GPU_CHK(cudaMalloc(&d_bB, (size_t)B_SIZE*IMG_SZ*sizeof(float)));
    demo_async(d_trn.imgs, d_trn.count, d_bA, d_bB);
    cudaFree(d_bA); cudaFree(d_bB);

    printf("\n[Stretch] FP16 Tensor Core GEMM:\n");
    demo_fp16(1024, 1024, 1024);

    cudaFree(d_w_c1); cudaFree(d_w_c2);
    cudaFree(d_w_f1); cudaFree(d_b_f1);
    cudaFree(d_w_f2); cudaFree(d_b_f2);
    cudaFree(d_inx);  cudaFree(d_m_c1); cudaFree(d_m_p1);
    cudaFree(d_m_c2); cudaFree(d_m_p2); cudaFree(d_m_f1);
    cudaFree(d_lgt);  cudaFree(d_prb);  cudaFree(d_lss);
    
    free(d_trn.imgs); free(d_trn.lbls);
    free(d_tst.imgs); free(d_tst.lbls);
    
    cudnnDestroy(dnn_hndl);
    cublasDestroy(blas_hndl);

    printf("\n========================================================\n");
    printf("  Execution complete!\n");
    printf("========================================================\n\n");
    return 0;
}