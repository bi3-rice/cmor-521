#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#define BLOCKSIZE 128

__global__ void stencil(const int N, const float alpha, float *y, const float *x){

  __shared__ float s_x[BLOCKSIZE+2];

  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int tid = threadIdx.x;
  s_x[tid] = 0.f;
  if (i < N){
    s_x[tid+1] = x[i];
  }
  if(tid == 0){
    if(i == 0){
      s_x[tid] = x[0];
    }
    else{
      s_x[tid] = x[i-1];
    }
  }
  else if(tid+1 == BLOCKSIZE){
    if(i == N-1){
      s_x[tid+2] == x[N-1];
    }
    else{
      s_x[tid+2] == x[i+1];
    }
  }
  
  __syncthreads();
  if (i < N){
    y[i] = alpha*(2*s_x[i+1] - x[i] - xM1[i+2]);
  }

}
    
int main(int argc, char * argv[]){

  int N = 4096;
  if (argc > 1){
    N = atoi(argv[1]);
  }

  int blockSize = BLOCKSIZE;
  float alpha = 1.f;

  // Next largest multiple of blockSize
  int numBlocks = (N + blockSize - 1) / blockSize;

  printf("N = %d, blockSize = %d, numBlocks = %d\n", N, blockSize, numBlocks);

  float * x = new float[N];
  float * y = new float[N];

  for (int i = 0; i < N; ++i){
    x[i] = 1.f;
  }

  // allocate memory and copy to the GPU
  float * d_x;
  float * d_y;  
  int size_xy = N * sizeof(float);
  int size_x_reduced = numBlocks * sizeof(float);
  cudaMalloc((void **) &d_x, size_xy);
  cudaMalloc((void **) &d_y, size_xy);
  
  // copy memory over to the GPU
  cudaMemcpy(d_x, x, size_xy, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, size_xy, cudaMemcpyHostToDevice);

  stencil <<< numBlocks, blockSize >>> (N, alpha, d_y, d_x);

  // copy memory back to the CPU
  cudaMemcpy(y, d_y, size_xy, cudaMemcpyDeviceToHost);

  bool equal = true;
  for (int i = 0; i < N; ++i){
    if(x[i] != 1 || y[i] != 0){
      equal = false;
      break;
    }
  }
  if(equal){
    printf("x[i] = 1, y[i] = 0 is true");
  }
  else{
    printf("x[i] = 1, y[i] = 0 is false");
  }

#if 1
  int num_trials = 10;
  float time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  for (int i = 0; i < num_trials; ++i){
    stencil <<< numBlocks, blockSize >>> (N, alpha, d_y, d_x);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  
  printf("Time to run kernel 10x: %6.2f ms.\n", time);
  
#endif

  return 0;
}
