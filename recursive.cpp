#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

#define BLOCK_SIZE 64

// computes C = C + A*B
void recur_2x2(double* C, double* A, double* B,
	       int CR_s, int CR_e, int CC_s, int CC_e,
	       int AR_s, int AR_e, int AC_s, int AC_e,
	       int BR_s, int BR_e, int BC_s, int BC_e, int n, int cur_n){
  if (cur_n < BLOCK_SIZE){
    int A_i = AR_s;
    int B_i = BR_s;
    for(int i = CR_s; i < CR_e; i++){
      int A_j = AC_s;
      int B_j = BC_s;
      for(int j = CC_s; j < CC_e; j++){
	C[i*n + j] = C[i*n + j] + A[A_i*n + A_j]*B[B_i*n + B_j];
	A_j += 1;
	B_j += 1;
      }
      A_i += 1;
      B_i += 1;
    }
    // C[CR_s*n + CC_s] = C[CR_s*n + CC_s] + A[AR_s*n + AC_s]*B[BR_s*n + BC_s];
  }
  else{
    
    // C00
    recur_2x2(C, A, B,
	      CR_s, (CR_e+CR_s)/2, CC_s, (CC_s+CC_e)/2, // C00
	      AR_s, (AR_s+AR_e)/2, AC_s, (AC_s+AC_e)/2, // A00
	      BR_s, (BR_s+BR_e)/2, BC_s, (BC_s+BC_e)/2, n, cur_n/2); // B00
    recur_2x2(C, A, B,
	      CR_s, (CR_s+CR_e)/2, CC_s, (CC_s+CC_e)/2, // C00
	      AR_s, (AR_s+AR_e)/2, (AC_e+AC_s)/2, AC_e, // A01
	      (BR_e+BR_s)/2, BR_e, BC_s, (BC_s+BC_e)/2, n, cur_n/2); // B10
    // C01
    recur_2x2(C, A, B,
	      CR_s, (CR_s+CR_e)/2, (CC_e+CC_s)/2, CC_e, // C01
	      AR_s, (AR_s+AR_e)/2, AC_s, (AC_s+AC_e)/2, // A00
	      BR_s, (BR_s+BR_e)/2, (BC_e+BC_s)/2, BC_e, n, cur_n/2); // B01
    recur_2x2(C, A, B,
	      CR_s, (CR_s+CR_e)/2, (CC_e+CC_s)/2, CC_e, // C01
	      AR_s, (AR_s+AR_e)/2, (AC_s+AC_e)/2, AC_e, // A01
	      (BR_s+BR_e)/2, BR_e, (BC_s+BC_e)/2, BC_e, n, cur_n/2); // B11
    // C10
    recur_2x2(C, A, B,
	      (CR_s+CR_e)/2, CR_e, CC_s, (CC_s+CC_e)/2, // C10
	      (AR_s+AR_e)/2, AR_e, AC_s, (AC_s+AC_e)/2, // A10
	      BR_s, (BR_s+BR_e)/2, BC_s, (BC_s+BC_e)/2, n, cur_n/2); // B00
    recur_2x2(C, A, B,
	      (CR_s+CR_e)/2, CR_e, CC_s, (CC_s+CC_e)/2, // C10
	      (AR_s+AR_e)/2, AR_e, (AC_s+AC_e)/2, AC_e, // A11
	      (BR_s+BR_e)/2, BR_e, BC_s, (BC_s+BC_e)/2, n, cur_n/2); // B10
    // C11
    recur_2x2(C, A, B,
	      (CR_s+CR_e)/2, CR_e, (CC_s+CC_e)/2, CC_e, // C11
	      (AR_s+AR_e)/2, AR_e, AC_s, (AC_s+AC_e)/2, // A10
	      BR_s, (BR_s+BR_e)/2, (BC_s+BC_e)/2, BC_e, n, cur_n/2); // B01
    recur_2x2(C, A, B,
	      (CR_s+CR_e)/2, CR_e, (CC_s+CC_e)/2, CC_e, // C11
	      (AR_s+AR_e)/2, AR_e, (AC_s+AC_e)/2, AC_e, // A11
	      (BR_s+BR_e)/2, BR_e, (BC_s+BC_e)/2, BC_e, n, cur_n/2); // B11
  }
}

/*
void strassen(double * C, double* A, double* B,
	      int CR_s, int CR_e, int CC_s, int CC_e,
	      int AR_s, int AR_e, int AC_s, int AC_e,
	      int BR_s, int BR_e, int BC_s, int BC_e, int n, int cur_n){
  if (cur_n == 2){
    cout << "Error? (if)" << endl;
    
    double M1 = (A[AR_s*n + AC_s] + A[AR_e*n + AC_e])*
      (B[BR_s*n + BC_s] + B[BR_e*n + BC_e]);
    
    double M2 = (A[AR_e*n + AC_s] + A[AR_e*n + AC_e])*
      B[BR_s*n + BC_s];
    
    double M3 = A[AR_s*n + AC_s]*
      (B[BR_s*n + BC_e] - B[BR_e*n + BC_e]);
    
    double M4 = A[AR_e*n + AC_e]*
      (B[BR_e*n + BC_s] - B[BR_s*n + BC_s]);
    
    double M5 = (A[AR_s*n + AC_s] + A[AR_s*n + AC_e])*
      B[BR_e*n + BC_e];
    
    double M6 = (A[AR_e*n + AC_s] - A[AR_s*n + AC_s])*
      (B[BR_s*n + BC_s] + B[BR_s*n + BC_e]);
    
    double M7 = (A[AR_s*n + AC_e] - A[AR_e*n + AC_e])*
      (B[BR_e*n + BC_s] + B[BR_e*n + BC_e]);
    
    /*
    double M1 = (A[0] + A[3])*(B[0] + B[3]);
    
    double M2 = (A[2] + A[3])*B[0];
    
    double M3 = A[0]*(B[1] - B[3]);
    
    double M4 = A[3]*(B[2] - B[0]);
    
    double M5 = (A[0] + A[1])*B[3];
    
    double M6 = (A[2] - A[0])*(B[0] + B[1]);
    
    double M7 = (A[1] - A[3])*(B[2] + B[3]);
    
    

    C[CR_s*n + CC_s] = M1 + M4 - M5 + M7;
    C[CR_s*n + CC_e] = M3 + M5;
    C[CR_e*n + CC_s] = M2 + M4;
    C[CR_e*n + CC_e] = M1 - M2 + M3 + M6;
    // C[CR_s*n + CC_s] = C[CR_s*n + CC_s] + A[AR_s*n + AC_s]*B[BR_s*n + BC_s];
  }
  else{
    cout << "Error? (else)" << endl;
    cout << "cur_n value:" << cur_n << endl;
    double * A11 = new double[cur_n * cur_n / 4];
    double * A12 = new double[cur_n * cur_n / 4];
    double * A21 = new double[cur_n * cur_n / 4];
    double * A22 = new double[cur_n * cur_n / 4];
    double * B11 = new double[cur_n * cur_n / 4];
    double * B12 = new double[cur_n * cur_n / 4];
    double * B21 = new double[cur_n * cur_n / 4];
    double * B22 = new double[cur_n * cur_n / 4];

    for(int i = 0; i < cur_n/2; i++){
      for(int j = 0; j < cur_n/2; j++){
	A11[i*n/2 + j] = A[i*cur_n/2 + j];
	A12[i*n/2 + j] = A[i*cur_n/2 + j + cur_n/2];
	A21[i*n/2 + j] = A[(i+1)*cur_n/2 + j];
	A22[i*n/2 + j] = A[(i+1)*cur_n/2 + j + cur_n/2];
	
	B11[i*n/2 + j] = B[i*cur_n/2 + j];
	B12[i*n/2 + j] = B[i*cur_n/2 + j + cur_n/2];
	B21[i*n/2 + j] = B[(i+1)*cur_n/2 + j];
	B22[i*n/2 + j] = B[(i+1)*cur_n/2 + j + cur_n/2];
      }
    }
    double * S1 = new double[cur_n * cur_n / 4];
    double * S2 = new double[cur_n * cur_n / 4];
    double * S3 = new double[cur_n * cur_n / 4];
    double * S4 = new double[cur_n * cur_n / 4];
    double * S5 = new double[cur_n * cur_n / 4];
    double * S6 = new double[cur_n * cur_n / 4];
    double * S7 = new double[cur_n * cur_n / 4];
    double * S8 = new double[cur_n * cur_n / 4];
    double * S9 = new double[cur_n * cur_n / 4];
    double * S10 = new double[cur_n * cur_n / 4];

    for(int i = 0; i < cur_n/2; i++){
      for(int j = 0; j < cur_n/2; j++){
        S1[cur_n*i + j] = A11[cur_n*i + j] + A22[cur_n*i + j];
	S2[cur_n*i + j] = B11[cur_n*i + j] + B22[cur_n*i + j];
	S3[cur_n*i + j] = A21[cur_n*i + j] + A22[cur_n*i + j];
	S4[cur_n*i + j] = B12[cur_n*i + j] - B22[cur_n*i + j];
	S5[cur_n*i + j] = B21[cur_n*i + j] - B11[cur_n*i + j];
	S6[cur_n*i + j] = A11[cur_n*i + j] + A12[cur_n*i + j];
	S7[cur_n*i + j] = A21[cur_n*i + j] - A11[cur_n*i + j];
	S8[cur_n*i + j] = B11[cur_n*i + j] + B12[cur_n*i + j];
	S9[cur_n*i + j] = A12[cur_n*i + j] - A22[cur_n*i + j];
	S10[cur_n*i + j] = B21[cur_n*i + j] + B22[cur_n*i + j];
      }
    }
    
    cout << "after for loop" << endl;
    // C00
    strassen(C, S1, S2,
			   CR_s, (CR_e+CR_s)/2, CC_s, (CC_s+CC_e)/2,
			   AR_s, (AR_s+AR_e)/2, AC_s, (AC_s+AC_e)/2,
			   BR_s, (BR_s+BR_e)/2, BC_s, (BC_s+BC_e)/2, n, cur_n/2);
    strassen(C, S3, B11,
			   CR_s, (CR_s+CR_e)/2, CC_s, (CC_s+CC_e)/2,
			   AR_s, (AR_s+AR_e)/2, (AC_e+AC_s)/2, AC_e,
			   (BR_e+BR_s)/2, BR_e, BC_s, (BC_s+BC_e)/2, n, cur_n/2);
    strassen(C, A11, S4,
			   CR_s, (CR_s+CR_e)/2, (CC_e+CC_s)/2, CC_e,
			   AR_s, (AR_s+AR_e)/2, AC_s, (AC_s+AC_e)/2,
			   BR_s, (BR_s+BR_e)/2, (BC_e+BC_s)/2, BC_e, n, cur_n/2);
    strassen(C, A22, S5,
			   CR_s, (CR_s+CR_e)/2, (CC_e+CC_s)/2, CC_e,
			   AR_s, (AR_s+AR_e)/2, (AC_s+AC_e)/2, AC_e,
			   (BR_s+BR_e)/2, BR_e, (BC_s+BC_e)/2, BC_e, n, cur_n/2);
    strassen(C, S6, B22,
			   (CR_s+CR_e)/2, CR_e, CC_s, (CC_s+CC_e)/2,
			   (AR_s+AR_e)/2, AR_e, AC_s, (AC_s+AC_e)/2,
			   BR_s, (BR_s+BR_e)/2, BC_s, (BC_s+BC_e)/2, n, cur_n/2); 
    strassen(C, S7, S8,
			   (CR_s+CR_e)/2, CR_e, CC_s, (CC_s+CC_e)/2,
			   (AR_s+AR_e)/2, AR_e, (AC_s+AC_e)/2, AC_e,
			   (BR_s+BR_e)/2, BR_e, BC_s, (BC_s+BC_e)/2, n, cur_n/2);
    strassen(C, S9, S10,
			  (CR_s+CR_e)/2, CR_e, (CC_s+CC_e)/2, CC_e,
			  (AR_s+AR_e)/2, AR_e, AC_s, (AC_s+AC_e)/2,
			   BR_s, (BR_s+BR_e)/2, (BC_s+BC_e)/2, BC_e, n, cur_n/2);
  }
}*/

double * strassen(double* A, double* B, int n, int cur_n){
  if (cur_n == 2){
    
    cout << "A Matrix:" << endl;
    for(int i = 0; i < cur_n; i++){
      cout << "[";
      for(int j = 0; j < cur_n; j++){
        cout << A[i*cur_n + j] << ", ";
      }
      cout << "]" << endl;
    }
    
    cout << "B Matrix:" << endl;
    for(int i = 0; i < cur_n; i++){
      cout << "[";
      for(int j = 0; j < cur_n; j++){
        cout << B[i*cur_n + j] << ", ";
      }
      cout << "]" << endl;
    }
    
    double M1 = (A[0] + A[3])*(B[0] + B[3]);
    
    double M2 = (A[2] + A[3])*B[0];
    
    double M3 = A[0]*(B[1] - B[3]);
    
    double M4 = A[3]*(B[2] - B[0]);
    
    double M5 = (A[0] + A[1])*B[3];
    
    double M6 = (A[2] - A[0])*(B[0] + B[1]);
    
    double M7 = (A[1] - A[3])*(B[2] + B[3]);
    
    cout << "right before declaring C" << endl;
    double * C = new double[cur_n * cur_n];
    cout << "right after declaring C" << endl;
      
    C[0] = M1 + M4 - M5 + M7;
    C[1] = M3 + M5;
    C[2] = M2 + M4;
    C[3] = M1 - M2 + M3 + M6;
    cout << "Right before return" << endl;
    return C;
    // C[CR_s*n + CC_s] = C[CR_s*n + CC_s] + A[AR_s*n + AC_s]*B[BR_s*n + BC_s];
  }
  else{
    cout << "Error? (else)" << endl;
    //cout << "cur_n value:" << cur_n << endl;
    double * A11 = new double[cur_n * cur_n / 4];
    double * A12 = new double[cur_n * cur_n / 4];
    double * A21 = new double[cur_n * cur_n / 4];
    double * A22 = new double[cur_n * cur_n / 4];
    double * B11 = new double[cur_n * cur_n / 4];
    double * B12 = new double[cur_n * cur_n / 4];
    double * B21 = new double[cur_n * cur_n / 4];
    double * B22 = new double[cur_n * cur_n / 4];

    for(int i = 0; i < cur_n/2; i++){
      for(int j = 0; j < cur_n/2; j++){
	A11[i*cur_n/2 + j] = A[i*cur_n/2 + j];
	A12[i*cur_n/2 + j] = A[i*cur_n/2 + j + cur_n/2];
	A21[i*cur_n/2 + j] = A[(i+1)*cur_n/2 + j];
	A22[i*cur_n/2 + j] = A[(i+1)*cur_n/2 + j + cur_n/2];
	
	B11[i*cur_n/2 + j] = B[i*cur_n/2 + j];
	B12[i*cur_n/2 + j] = B[i*cur_n/2 + j + cur_n/2];
	B21[i*cur_n/2 + j] = B[(i+1)*cur_n/2 + j];
	B22[i*cur_n/2 + j] = B[(i+1)*cur_n/2 + j + cur_n/2];
      }
    }
    double * S1 = new double[cur_n * cur_n / 4];
    double * S2 = new double[cur_n * cur_n / 4];
    double * S3 = new double[cur_n * cur_n / 4];
    double * S4 = new double[cur_n * cur_n / 4];
    double * S5 = new double[cur_n * cur_n / 4];
    double * S6 = new double[cur_n * cur_n / 4];
    double * S7 = new double[cur_n * cur_n / 4];
    double * S8 = new double[cur_n * cur_n / 4];
    double * S9 = new double[cur_n * cur_n / 4];
    double * S10 = new double[cur_n * cur_n / 4];

    for(int i = 0; i < cur_n/2; i++){
      for(int j = 0; j < cur_n/2; j++){
        S1[cur_n*i + j] = A11[cur_n*i + j] + A22[cur_n*i + j];
	S2[cur_n*i + j] = B11[cur_n*i + j] + B22[cur_n*i + j];
	S3[cur_n*i + j] = A21[cur_n*i + j] + A22[cur_n*i + j];
	S4[cur_n*i + j] = B12[cur_n*i + j] - B22[cur_n*i + j];
	S5[cur_n*i + j] = B21[cur_n*i + j] - B11[cur_n*i + j];
	S6[cur_n*i + j] = A11[cur_n*i + j] + A12[cur_n*i + j];
	S7[cur_n*i + j] = A21[cur_n*i + j] - A11[cur_n*i + j];
	S8[cur_n*i + j] = B11[cur_n*i + j] + B12[cur_n*i + j];
	S9[cur_n*i + j] = A12[cur_n*i + j] - A22[cur_n*i + j];
	S10[cur_n*i + j] = B21[cur_n*i + j] + B22[cur_n*i + j];
      }
    }

    cout << "first strassen" << endl;
    double * M1 = strassen(S1, S2, n, cur_n/2);
    cout << "second strassen" << endl;
    double * M2 = strassen(S3, B11, n, cur_n/2);
    cout << "third  strassen" << endl;
    double * M3 = strassen(A11, S4, n, cur_n/2);
    cout << "fourth strassen" << endl;
    double * M4 = strassen(A22, S5, n, cur_n/2);
    cout << "fifth strassen" << endl;
    double * M5 = strassen(S6, B22, n, cur_n/2);
    cout << "sixth strassen" << endl;
    double * M6 = strassen(S7, S8, n, cur_n/2);
    cout << "seventh strassen" << endl;
    double * M7 = strassen(S9, S10, n, cur_n/2);

    double * C = new double[cur_n * cur_n];

    /*
    for(int i = 0; i < cur_n/2; i++){
      for(int j = 0; j < cur_n/2; j++){
	// C00
	C[i*cur_n/2 + j] = M1[i*cur_n/2 + j] + M4[i*cur_n/2 + j] - M5[i*cur_n/2 + j] + M7[i*cur_n/2 + j];
	// C01
	C[i*cur_n/2 + j + cur_n/2] = M3[i*cur_n/2 + j] + M5[i*cur_n/2 + j];
	// C10
	C[(i+cur_n/2)*cur_n/2 + j] = M2[i*cur_n/2 + j] + M4[i*cur_n/2 + j];
	// C11
	C[(i+cur_n/2)*cur_n/2 + j + cur_n/2] = M1[i*cur_n/2 + j] - M2[i*cur_n/2 + j] + M3[i*cur_n/2 + j] + M6[i*cur_n/2 + j];
      }
    }*/
    
    return C;
  }

  
}

void matmul_blocked(const int n, double* C, double* A, double* B){
  for (int i = 0; i < n; i += BLOCK_SIZE){
    for (int j = 0; j < n; j += BLOCK_SIZE){
      for (int k = 0; k < n; k += BLOCK_SIZE){

	// small matmul
	for (int ii = i; ii < i + BLOCK_SIZE; ii++){
	  for (int jj = j; jj < j + BLOCK_SIZE; jj++){
	    double Cij = C[jj + ii * n];
	    for (int kk = k; kk < k + BLOCK_SIZE; kk++){
	      Cij += A[kk + ii * n] * B[jj + kk * n]; // Aik * Bkj
	    }
	    C[jj + ii * n] = Cij;
	  }
	}
	
      }
    }
  }
}

int main(int argc, char * argv[]){

  int n = atoi(argv[1]);
  cout << "Matrix size n = " << n << ", block size = " << BLOCK_SIZE << endl;
  
  double * A = new double[n * n];
  double * B = new double[n * n];
  double * C = new double[n * n];

  // make A, B = I
  for (int i = 0; i < n; ++i){
    A[i + i * n] = 1.0;
    B[i + i * n] = 1.0;
  }
  

  int num_trials = 5;

  // Measure performance
  high_resolution_clock::time_point start = high_resolution_clock::now();
  for (int i = 0; i < num_trials; ++i){
    for (int i = 0; i < n * n; ++i){
      C[i] = 0.0;
    }
    recur_2x2(C, A, B,
	      0, n, 0, n,
	      0, n, 0, n,
	      0, n, 0, n, n, n);
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  duration<double> elapsed_naive = (end - start) / num_trials;

  double sum_C = 0.0;
  for (int i = 0; i < n * n; ++i){
    sum_C += C[i];
  }
  cout << "Recur_2x2 sum_C = " << sum_C << endl;

  cout << "Recur_2x2 elapsed time (ms) = " << elapsed_naive.count() * 1000 << endl;

  start = high_resolution_clock::now();
  for (int i = 0; i < num_trials; ++i){
    for (int i = 0; i < n * n; ++i){
      C[i] = 0.0;
    }
    C = strassen(A, B, n, n);
  }
  end = high_resolution_clock::now();
  duration<double> elapsed_strass = (end - start) / num_trials;

  sum_C = 0.0;
  for (int i = 0; i < n * n; ++i){
    sum_C += C[i];
  }
  cout << "Strassen_2x2 sum_C = " << sum_C << endl;

  cout << "Strassen elapsed time (ms) = " << elapsed_strass.count() * 1000 << endl;

  /*
  cout << "Matrix C:" << endl;
  
  for(int i = 0; i < n; i++){
    cout << "[";
    for(int j = 0; j < n; j++){
      cout << C[i*n + j] << ", ";
    }
    cout << "]" << endl;
  }

  cout << "Matrix A:" << endl;
  
  for(int i = 0; i < n; i++){
    cout << "[";
    for(int j = 0; j < n; j++){
      cout << A[i*n + j] << ", ";
    }
    cout << "]" << endl;
  }

  cout << "Matrix B:" << endl;
  
  for(int i = 0; i < n; i++){
    cout << "[";
    for(int j = 0; j < n; j++){
      cout << B[i*n + j] << ", ";
    }
    cout << "]" << endl;
  }
  */

  delete[] A;
  delete[] B;
  delete[] C;  
  
  return 0;
}
