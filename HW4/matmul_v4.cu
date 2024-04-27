#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#define BLOCKSIZE 32

__global__ void sgemm1DBlocktiling(int N, float alpha,
                                   const float *A, const float *B, float beta, float *C) {
  // If we flip x and y here we get ~30% less performance for large matrices.
  // The current, 30% faster configuration ensures that blocks with sequential
  // blockIDs access columns of B sequentially, while sharing the same row of A.
  // The slower configuration would share columns of A, but access into B would
  // be non-sequential. So the faster configuration has better spatial locality
  // and hence a greater L2 hit rate.
  const uint cRow = blockIdx.y;
  const uint cCol = blockIdx.x;
  const uint TM = 16;

  // each warp will calculate 32*TM elements, with 32 being the columnar dim.
  const int threadCol = threadIdx.x % BLOCKSIZE;
  const int threadRow = threadIdx.x / BLOCKSIZE;

  // allocate space for the current blocktile in SMEM
  __shared__ float As[BLOCKSIZE * BLOCKSIZE];
  __shared__ float Bs[BLOCKSIZE * BLOCKSIZE];

  // Move blocktile to beginning of A's row and B's column
  A += cRow * BLOCKSIZE * N;
  B += cCol * BLOCKSIZE;
  C += cRow * BLOCKSIZE * N + cCol * BLOCKSIZE;

  // todo: adjust this to each thread to load multiple entries and
  // better exploit the cache sizes
  // assert(BLOCKSIZE * BLOCKSIZE == blockDim.x);
  // assert(BLOCKSIZE * BLOCKSIZE == blockDim.x);
  const uint innerColA = threadIdx.x % BLOCKSIZE; // warp-level GMEM coalescing
  const uint innerRowA = threadIdx.x / BLOCKSIZE;
  const uint innerColB = threadIdx.x % BLOCKSIZE; // warp-level GMEM coalescing
  const uint innerRowB = threadIdx.x / BLOCKSIZE;

  // allocate thread-local cache for results in registerfile
  float * threadResults = new float[TM];

  // outer loop over block tiles
  for (uint bkIdx = 0; bkIdx < N; bkIdx += BLOCKSIZE) {
    // populate the SMEM caches
    As[innerRowA * BLOCKSIZE + innerColA] = A[innerRowA * N + innerColA];
    Bs[innerRowB * BLOCKSIZE + innerColB] = B[innerRowB * N + innerColB];
    __syncthreads();

    // advance blocktile
    A += BLOCKSIZE;
    B += BLOCKSIZE * N;

    // calculate per-thread results
    for (uint dotIdx = 0; dotIdx < BLOCKSIZE; ++dotIdx) {
      // we make the dotproduct loop the outside loop, which facilitates
      // reuse of the Bs entry, which we can cache in a tmp var.
      float tmpB = Bs[dotIdx * BLOCKSIZE + threadCol];
      for (uint resIdx = 0; resIdx < TM; ++resIdx) {
        threadResults[resIdx] +=
            As[(threadRow * TM + resIdx) * BLOCKSIZE + dotIdx] * tmpB;
      }
    }
    __syncthreads();
  }

  // write out the results
  for (uint resIdx = 0; resIdx < TM; ++resIdx) {
    C[(threadRow * TM + resIdx) * N + threadCol] =
        alpha * threadResults[resIdx] +
        beta * C[(threadRow * TM + resIdx) * N + threadCol];
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
  sgemm1DBlocktiling <<< gridDims, blockDims >>> (N, alpha, d_A, d_B, beta, d_C);

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
