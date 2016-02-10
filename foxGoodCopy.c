void Parallel_matrix_mult(float local_A[], float local_B[]),
float local_C[], int n,  int n_bar,  int p){
    
    float *        B_cols;
    MPI_Datatype  gather_mpi_t;
    int           block;
    
    Allocate matrix(&B_cols,  n,  n_bar,   "B_cols");
    MPI_Type_vector(n_bar,  n_bar,  n,  MPI_FLOAT, &gather_mpi_t);
    MPI_Type_commit(&gather_mpi_t);
    
    for  (block = 0; block < p;  block++){
        MPI_Allgather(local_B + block*n_bar,1, gather_mpi_t,
                      B_cols,  n_bar*n_bar,  MPI_FLOAT,  MPI_COMM_WORLD);
        Matrix_mult(local_A,  B_cols,  local_C,  n_bar, n, block);
    }
    free(B_cols);
    MPI_Type_free(&gather_mpi_t);
} /*Parallel matrix mult*/