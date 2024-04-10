#include "mpi.h"
#include <iostream>
#include <cmath>

using namespace std;

void matmul(double * A, double * B, double * C, int n){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      for(int k = 0; k < n; k++){
	C[j + i * n] += A[k + i * n] * B[j + k * n];
      }
    }
  }
}

bool check_equal(double * C1, double * C2, int n){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(abs(C1[j + i * n] - C2[j + i * n]) > 1e-6){
	return false;
	break;
      }
    }
  }
  return true;
}

int main(int argc, char* argv[]){
  int n = 128;
  if(argc > 0){
    n = atoi(argv[1]);
  }
  int p;

  double * A = new double[n * n];
  double * B = new double[n * n];
  double * C = new double[n * n];
  
  MPI_Init(NULL, NULL);
  MPI_Status status;
  int rank, num_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  for(int i = 0; i < n * n; i++){
    A[i] = rand() % 4;
    B[i] = rand() % 4;
    C[i] = 0.0;
  }

  p = (int) sqrt(num_procs);
  
  matmul(A, B, C, n);

  if(rank == 0 && n <= 16){
    cout << "Matrix A" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << A[j + i*n] << ",";
      }
      cout << "]" << endl;
    }
    
    cout << "Matrix B" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << B[j + i*n] << ",";
      }
      cout << "]" << endl;
    }

    cout << "Matrix C" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << C[j + i*n] << ",";
      }
      cout << "]" << endl;
    }
  }
  //MPI_Barrier(MPI_COMM_WORLD);
  
  int sub_n = n / p;

  int row = rank / p;
  int col = rank % p;
  int a_start_row = (row)*sub_n;
  int k = ((row + col) % p)*sub_n;
  int b_start_col = (col)*sub_n;

  int c_start_row = a_start_row;
  int c_start_col = b_start_col;
  
  double * first_sub_A = new double[sub_n * sub_n];
  double * first_sub_B = new double[sub_n * sub_n];
  
  double * sub_C = new double[sub_n * sub_n];
  double * rec_C = new double[sub_n * sub_n];
  int rec_row = 0;
  int rec_col = 0;


  for(int i = 0; i < sub_n; i++){
    for(int j = 0; j < sub_n; j++){
      first_sub_A[j + i*sub_n] = A[(j + k) + (i + a_start_row)*n];
      first_sub_B[j + i*sub_n] = B[(j + b_start_col) + (i + k)*n];
      sub_C[j + i*sub_n] = 0.0;
    }
  }

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, row, rank, &row_comm);

  MPI_Comm col_comm;
  MPI_Comm_split(MPI_COMM_WORLD, col, rank, &col_comm);

  matmul(first_sub_A, first_sub_B, sub_C, sub_n);
  int left_rank, right_rank, above_rank, below_rank;

  left_rank = (col - 1 + p) % p;
  right_rank = (col + 1) % p;

  above_rank = (row - 1 + p) % p;
  below_rank = (row + 1) % p;

  
  for(int i = 0; i < p-1; i++){    
    
    MPI_Sendrecv_replace(first_sub_A, sub_n*sub_n, MPI_DOUBLE, left_rank, 1, right_rank, 1, row_comm, &status);
    MPI_Sendrecv_replace(first_sub_B, sub_n*sub_n, MPI_DOUBLE, above_rank, 1, below_rank, 1, col_comm, &status);

    MPI_Barrier(MPI_COMM_WORLD);
    
    matmul(first_sub_A, first_sub_B, sub_C, sub_n);
    
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  double * calc_C = new double[n * n];

  if(rank == 0){
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	calc_C[j + i * n] = sub_C[j + i * sub_n];
      }
    }
  }

  for(int id = 1; id < num_procs; id++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == id){
      MPI_Send(sub_C, sub_n * sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&c_start_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&c_start_col, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
    if(rank == 0){
      MPI_Recv(rec_C, sub_n * sub_n, MPI_DOUBLE, id, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&rec_row, 1, MPI_INT, id, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&rec_col, 1, MPI_INT, id, 1, MPI_COMM_WORLD, &status);
      for(int i = 0; i < sub_n; i++){
	for(int j = 0; j < sub_n; j++){
	  calc_C[(j + rec_col) + (i + rec_row)*n] = rec_C[j + i * sub_n];
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(rank == 0){
    if(n <= 16){
      cout << "Calculated C, n = " << n << endl;
      for(int i = 0; i < n; i++){
        cout << "[";
        for(int j = 0; j < n; j++){
	  cout << calc_C[j + i * n] << ",";
        }
        cout << "]" << endl;
      }
    }

    bool equal = check_equal(C, calc_C, n);
    if(equal){
      cout << "C and Cannon's C are equal" << endl;
    }
    else{
      cout << "C and Cannon's C are not equal" << endl;
    }
  }
  
  delete[] first_sub_A;
  delete[] first_sub_B;
  delete[] sub_C;
  delete[] rec_C;
  delete[] calc_C;
  
  MPI_Finalize();

  delete[] C;
  delete[] A;
  delete[] B;

  return 0;
}
