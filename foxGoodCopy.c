void Parallel_matrix_mult(float local_A[], float local_B[]),
float local_C[], int n,  int n_bar,  int p){
    
    float*        B_cols;
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

void Setup_grid(MPI_Comm*  grid_comm,  MPI_Comm*  row_comm,
                MPI_comm*  col_comm,  int*  my_row,  int*  my_col,
                int*  my_grid_rank){
    int old_rank,  dimensions[2],     wrap_around[2]
    int coordinates[2],     free_coords[2],     p,  q,  old_rank;

    MPI_Comm_size(MPI_COMM_WORLD,   &p);
    MPI_Comm_rank(MPI_COMM_WORLD,   &old_rank);
    
    q = (int)  sqrt((double)  p);
    dimensions[0] =  dimensions[1] =  q;
    
wrap_around[0]      wrap around
MPI Cart create MPI COMM WORLD     dimensions
wrap around      grid comm
MPI Comm rank grid comm   my grid rank
MPI Cart coords grid comm  my grid rank
coordinates
my row   coordinates
my col   coordinates
free coords         free coords
MPI Cart sub grid comm  free coords   row comm
free coords         free coords
MPI Cart sub grid comm  free coords   col comm
Setup grid
}