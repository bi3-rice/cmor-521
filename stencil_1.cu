#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

#define BLOCKSIZE 256

__global__ void stencil(const int N, const float alpha, float *y, const float *x){

  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int tid = threadIdx.x;
  
  if (i < N){
    float xM1;
    float xP1;
    
    if (i == 0){
      xM1 = x[0];
      xP1 = x[i+1];
    }
    else if(i == N-1){
      xM1 = x[i-1];
      xP1 = x[N-1];
    }
    else{
      xM1 = x[i-1];
      xP1 = x[i+1];
    }
    y[i] = alpha*(2*x[i] - xP1 - xM1);
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
