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
  
  int n = atoi(argv[1]);
  int p = 2;

  double * A = new double[n * n];
  double * B = new double[n * n];
  double * C = new double[n * n];
  
  
  // since p = 2, the below code assumes 4 processors
  
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

  matmul(A, B, C, n);

  if(rank == 0){
    cout << "Matrix A" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << A[j + i*n] << ",";
      }
      cout << "]" << endl;
    }
    cout << endl;
    
    cout << "Matrix B" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << B[j + i*n] << ",";
      }
      cout << "]" << endl;
    }
    cout << endl;

    cout << "Matrix C" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << C[j + i*n] << ",";
      }
      cout << "]" << endl;
    }
    cout << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  int sub_n = n / p;

  int row = rank / p;
  int col = rank % p;
  int start_row = (row)*sub_n;
  int start_col = (col)*sub_n;
  
  double * first_sub_A = new double[sub_n * sub_n];
  double * first_sub_B = new double[sub_n * sub_n];
  double * second_sub_A = new double[sub_n * sub_n];
  double * second_sub_B = new double[sub_n * sub_n];
  
  bool first_A_def, first_B_def;
  double * sub_C = new double[sub_n * sub_n];
  double * rec_C = new double[sub_n * sub_n];
  int rec_row = 0;
  int rec_col = 0;

  if(rank == 0){
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	first_sub_A[j + i*sub_n] = A[(j + start_col) + (i + start_row)*n];
	second_sub_A[j + i*sub_n] = 0.0;
	first_sub_B[j + i*sub_n] = B[(j + start_col) + (i + start_row)*n];
	second_sub_B[j + i*sub_n] = 0.0;
	sub_C[j + i*sub_n] = 0.0;
      }
    }
    first_A_def = true;
    first_B_def = true;
  }
  if(rank == 1){
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	first_sub_A[j + i*sub_n] = 0.0;
	second_sub_A[j + i*sub_n] = A[(j + start_col) + (i + start_row)*n];
	first_sub_B[j + i*sub_n] = B[(j + start_col) + (i + start_row)*n];
	second_sub_B[j + i*sub_n] = 0.0;
	sub_C[j + i*sub_n] = 0.0;
      }
    }
    first_A_def = false;
    first_B_def = true;
  }
  if(rank == 2){
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	first_sub_A[j + i*sub_n] = A[(j + start_col) + (i + start_row)*n];
	second_sub_A[j + i*sub_n] = 0.0;
	first_sub_B[j + i*sub_n] = 0.0;
	second_sub_B[j + i*sub_n] = B[(j + start_col) + (i + start_row)*n];
	sub_C[j + i*sub_n] = 0.0;
      }
    }
    first_A_def = true;
    first_B_def = false;
  }
  if(rank == 3){
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	first_sub_A[j + i*sub_n] = 0.0;
	second_sub_A[j + i*sub_n] = A[(j + start_col) + (i + start_row)*n];
	first_sub_B[j + i*sub_n] = 0.0;
	second_sub_B[j + i*sub_n] = B[(j + start_col) + (i + start_row)*n];
	sub_C[j + i*sub_n] = 0.0;
      }
    }
    first_A_def = false;
    first_B_def = false;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, row, rank, &row_comm);

  MPI_Comm col_comm;
  MPI_Comm_split(MPI_COMM_WORLD, col, rank, &col_comm);

  for(int i = 0; i < p; i++){
    if(i == 0){
      MPI_Bcast(first_sub_A, sub_n*sub_n, MPI_DOUBLE, i, row_comm);
      MPI_Bcast(first_sub_B, sub_n*sub_n, MPI_DOUBLE, i, col_comm);
    }
    else{
      MPI_Bcast(second_sub_A, sub_n*sub_n, MPI_DOUBLE, i, row_comm);
      MPI_Bcast(second_sub_B, sub_n*sub_n, MPI_DOUBLE, i, col_comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  matmul(first_sub_A, first_sub_B, sub_C, sub_n);
  matmul(second_sub_A, second_sub_B, sub_C, sub_n);

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
      MPI_Send(&start_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&start_col, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
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
    cout << "Calculated C" << endl;
	
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << calc_C[j + i * n] << ",";
      }
      cout << "]" << endl;
    }

    bool equal = check_equal(C, calc_C, n);
    if(equal){
      cout << "C and SUMMA C are equal" << endl;
    }
    else{
      cout << "C and SUMMA C are not equal" << endl;
    }
  }
  
  delete[] first_sub_A;
  delete[] first_sub_B;
  delete[] second_sub_A;
  delete[] second_sub_B;
  delete[] sub_C;
  delete[] rec_C;
  delete[] calc_C;
  
  MPI_Finalize();

  delete[] C;
  delete[] A;
  delete[] B;

  return 0;
}
