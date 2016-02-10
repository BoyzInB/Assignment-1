#include "matrix.h"

int main()
{
    matrix_t matrixA, matrixB, matrix;
    int n;
//    double max_norm;
    printf("Perform matrix-matrix multiply\n");

    // Prompt user and read matrix size from stdin.
    printf("Enter the Dimension for a square matrix: ");
    scanf("%d",&n);
    printf("allocated matrix %d * %d\n", n,n);

    matrix_allocate(&matrix, n);
    matrix_allocate(&matrixA,n);
    matrix_allocate(&matrixB,n);
    // Fill matrix.
    matrix_fill_random(matrixA);
    matrix_fill_random(matrixB);
    matrix_fill_zeros(matrix);
    // Compute max norm of matrix and print it.
    //max_norm = matrix_max_norm(matrix);
    matrix = matrix_mul(matrixA, matrixB,n);
    //matrix_deallocate(&matrix);
    //printf("Maximum norm : %.4f\n", max_norm);
    
    // Print matrix.
    printf("A is \n");
    matrix_print(matrixA);
    printf("B is \n");
    matrix_print(matrixB);
    printf("C = A*B  is \n");
    matrix_print(matrix);
    
    // Deallocate matrix.    
    matrix_deallocate(&matrix);
    matrix_deallocate(&matrixA);
    matrix_deallocate(&matrixB);

    // Exit.
    return 0;
}


