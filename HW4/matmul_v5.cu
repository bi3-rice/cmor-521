#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#define BLOCKSIZE 32

__global__ void sgemm2DBlocktiling(int N, float alpha, const float *A,
                       const float *B, float beta, float *C) {
  const uint cRow = blockIdx.y;
  const uint cCol = blockIdx.x;
  const uint TM = 4;
  const uint TN = 4;

  const uint totalResultsBlocktile = BLOCKSIZE * BLOCKSIZE;
  // A thread is responsible for calculating TM*TN elements in the blocktile
  const uint numThreadsBlocktile = totalResultsBlocktile / (TM * TN);

  // ResultsPerBlock / ResultsPerThread == ThreadsPerBlock
  // assert(numThreadsBlocktile == blockDim.x);

  // BN/TN are the number of threads to span a column
  const int threadCol = threadIdx.x % (BLOCKSIZE / TN);
  const int threadRow = threadIdx.x / (BLOCKSIZE / TN);

  // allocate space for the current blocktile in smem
  __shared__ float As[BLOCKSIZE * BLOCKSIZE];
  __shared__ float Bs[BLOCKSIZE * BLOCKSIZE];

  // Move blocktile to beginning of A's row and B's column
  A += cRow * BLOCKSIZE * N;
  B += cCol * BLOCKSIZE;
  C += cRow * BLOCKSIZE * N + cCol * BLOCKSIZE;

  // calculating the indices that this thread will load into SMEM
  const uint innerRowA = threadIdx.x / BLOCKSIZE;
  const uint innerColA = threadIdx.x % BLOCKSIZE;
  // calculates the number of rows of As that are being loaded in a single step
  // by a single block
  const uint strideA = numThreadsBlocktile / BLOCKSIZE;
  const uint innerRowB = threadIdx.x / BLOCKSIZE;
  const uint innerColB = threadIdx.x % BLOCKSIZE;
  // for both As and Bs we want each load to span the full column-width, for
  // better GMEM coalescing (as opposed to spanning full row-width and iterating
  // across columns)
  const uint strideB = numThreadsBlocktile / BLOCKSIZE;

  // allocate thread-local cache for results in registerfile
  float threadResults[TM * TN] = {0.0};
  // register caches for As and Bs
  float regM[TM] = {0.0};
  float regN[TN] = {0.0};

  // outer-most loop over block tiles
  for (uint bkIdx = 0; bkIdx < N; bkIdx += BLOCKSIZE) {
    // populate the SMEM caches
    for (uint loadOffset = 0; loadOffset < BLOCKSIZE; loadOffset += strideA) {
      As[(innerRowA + loadOffset) * BLOCKSIZE + innerColA] =
          A[(innerRowA + loadOffset) * N + innerColA];
    }
    for (uint loadOffset = 0; loadOffset < BLOCKSIZE; loadOffset += strideB) {
      Bs[(innerRowB + loadOffset) * BLOCKSIZE + innerColB] =
          B[(innerRowB + loadOffset) * N + innerColB];
    }
    __syncthreads();

    // advance blocktile
    A += BLOCKSIZE;     // move BK columns to right
    B += BLOCKSIZE * N; // move BK rows down

    // calculate per-thread results
    for (uint dotIdx = 0; dotIdx < BLOCKSIZE; ++dotIdx) {
      // block into registers
      for (uint i = 0; i < TM; ++i) {
        regM[i] = As[(threadRow * TM + i) * BLOCKSIZE + dotIdx];
      }
      for (uint i = 0; i < TN; ++i) {
        regN[i] = Bs[dotIdx * BLOCKSIZE + threadCol * TN + i];
      }
      for (uint resIdxM = 0; resIdxM < TM; ++resIdxM) {
        for (uint resIdxN = 0; resIdxN < TN; ++resIdxN) {
          threadResults[resIdxM * TN + resIdxN] +=
              regM[resIdxM] * regN[resIdxN];
        }
      }
    }
    __syncthreads();
  }

  // write out the results
  for (uint resIdxM = 0; resIdxM < TM; ++resIdxM) {
    for (uint resIdxN = 0; resIdxN < TN; ++resIdxN) {
      C[(threadRow * TM + resIdxM) * N + threadCol * TN + resIdxN] =
          alpha * threadResults[resIdxM * TN + resIdxN] +
          beta * C[(threadRow * TM + resIdxM) * N + threadCol * TN + resIdxN];
    }
  }
}


int main(int argc, char * argv[]){

  float alpha = 0.5, beta = 3.0;
  
  int N = 4096;
  if (argc > 1){
    N = atoi(argv[1]);
  }

  float * A = new float[N * N];
  float * B = new float[N * N];
  float * C = new float[N * N];

  for (int i = 0; i < N * N; ++i){
    A[i] = 0.f;
    B[i] = 0.f;
    C[i] = 0.f;
  }
  for (int i = 0; i < N; ++i){
    A[i + i * N] = 1.f; // identity
    B[i + i * N] = 1.f; // identity
  }

  // allocate memory and copy to the GPU
  float * d_A;
  float * d_B;
  float * d_C;
  int size = N * N * sizeof(float);
  cudaMalloc((void **) &d_A, size);
  cudaMalloc((void **) &d_B, size);
  cudaMalloc((void **) &d_C, size);
  
  // copy memory over to the GPU
  cudaMemcpy(d_A, A, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_C, C, size, cudaMemcpyHostToDevice);

  // Next largest multiple of blockSize
  int numBlocks = (N + BLOCKSIZE - 1) / BLOCKSIZE; 
  printf("N = %d, numBlocks * blockSize = %d\n", N, numBlocks * BLOCKSIZE);
  dim3 gridDims(numBlocks, numBlocks);
  dim3 blockDims(BLOCKSIZE * BLOCKSIZE);
  sgemm2DBlocktiling <<< gridDims, blockDims >>> (N, alpha, d_A, d_B, beta, d_C);

  // copy memory back to the CPU
  cudaMemcpy(C, d_C, size, cudaMemcpyDeviceToHost);
  
  float error = 0.f;
  for (int i = 0; i < N; ++i){
    for (int j = 0; j < N; ++j){
      //      printf("C[%d,%d] = %f\n", i, j, C[j + i * N]);
      float Cij = 0.f;
      if (i==j){
	Cij = 1.f;
      }
      float diff = C[j + i * N] - Cij;
      error += fabs(diff);
    }
  }
  printf("error = %f\n", error);


  return 0;
}
