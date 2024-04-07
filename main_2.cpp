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

int main(){
  
  int n = 4;
  int p = 2;

  double * A = new double[n * n];
  double * B = new double[n * n];
  
  // since p = 2, the below code assumes 4 processors
  
  MPI_Init(NULL, NULL);
  MPI_Status status;
  int rank, num_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  for(int i = 0; i < n * n; i++){
    A[i] = rand() % 4;
    B[i] = rand() % 4;
  }

  if(rank == 0){
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
    for(int j = 0; j < p; j++){
      if(i == row && j == col){
	if(first_A_def){
	  MPI_Bcast(first_sub_A, sub_n*sub_n, MPI_DOUBLE, i, row_comm);
	}
	else{
	  MPI_Bcast(second_sub_A, sub_n*sub_n, MPI_DOUBLE, i, row_comm);
	}

	if(first_B_def){
	  MPI_Bcast(first_sub_B, sub_n*sub_n, MPI_DOUBLE, j, col_comm);
	}
	else{
	  MPI_Bcast(second_sub_B, sub_n*sub_n, MPI_DOUBLE, j, col_comm);
	}
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  matmul(first_sub_A, first_sub_B, sub_C, sub_n);
  matmul(second_sub_A, second_sub_B, sub_C, sub_n);
  
  cout << endl;

  for(int id = 0; id < p*p; id++){
    if(rank == id){
      cout << "=================" << endl;
      cout << "===PROCESSOR " << rank << "===" << endl;
      cout << "=================" << endl;
      cout << "1st Submatrix A" << endl;
      
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << first_sub_A[j + i*sub_n] << ",";
	}
	cout << "]" << endl;
      }
      
      cout << "1st Submatrix B" << endl;
      
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << first_sub_B[j + i*sub_n] << ",";
	}
	cout << "]" << endl;
      }
      
      cout << "2nd Submatrix A" << endl;
      
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << second_sub_A[j + i*sub_n] << ",";
	}
	cout << "]" << endl;
      }
      
      cout << "2nd Submatrix B" << endl;
      
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << second_sub_B[j + i*sub_n] << ",";
	}
	cout << "]" << endl;
      }
      
      
      cout << "Submatrix C" << endl;
	
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << sub_C[j + i*sub_n] << ",";
	}
	cout << "]" << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Request send_req;

  /*
  if(rank == 1){
    MPI_Isend(sub_C12, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);
  }
  if(rank == 2){
    MPI_Isend(sub_C21, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);
  }
  if(rank == 3){
    MPI_Isend(sub_C22, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);
  }
  if(rank == 0){
    MPI_Irecv(sub_C12, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);
    MPI_Irecv(sub_C21, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);
    MPI_Irecv(sub_C22, sub_n*sub_n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &send_req);

    double * C = new double[n*n];
    
    for(int i = 0; i < sub_n; i++){
      for(int j = 0; j < sub_n; j++){
	C[j + i * n] = sub_C11[j + i * sub_n];
	C[(j+sub_n) + i * n] = sub_C12[j + i * sub_n];
	C[j + (i+sub_n) * n] = sub_C21[j + i * sub_n];
	C[(j+sub_n) + (i+sub_n) * n] = sub_C22[j + i * sub_n];
      }
    }

    cout << "Matrix C" << endl;
    for(int i = 0; i < n; i++){
      cout << "[";
      for(int j = 0; j < n; j++){
	cout << C[j + i * n] << ", ";
      }
      cout << "]" << endl;
    }

    delete[] C;
  }
  */

  /*
  for(int id = 0; id < p*p; id++){
    if(rank == id){
      cout << "Rank = " << id << endl;
      cout << "Submatrix C" << endl;
      
      for(int i = 0; i < sub_n; i++){
	cout << "[";
	for(int j = 0; j < sub_n; j++){
	  cout << sub_C[j + i * n] << ", ";
	}
	cout << "]" << endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  */
  
  delete[] first_sub_A;
  delete[] first_sub_B;
  delete[] second_sub_A;
  delete[] second_sub_B;
  
  MPI_Finalize();

  delete[] A;
  delete[] B;

  return 0;
}
