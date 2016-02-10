#include "matrix.h"

/**
 *  Allocate matrix data to given square size. 
 */
void matrix_allocate(matrix_t *matrix, int n)
{
    matrix->n = n;
    matrix->data = (double **) malloc(n*sizeof(double *));
    if(matrix->data == NULL) {
        printf("Unable to allocate matrix.");
        exit(-1);
    }
    int i;
    for(i = 0; i < n; i++) {
        matrix->data[i] = (double *) malloc(n*sizeof(double));
        if(matrix->data[i] == NULL) {
            printf("Unable to allocate matrix.");
            exit(-1);
        }
    }
}

/**
 *  Deallocate matrix data. 
 */
void matrix_deallocate(matrix_t *matrix)
{
    int i;
    for(i = 0; i < matrix->n; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
}

/** 
 *  Fill the matrix with random values.
 */
void matrix_fill_random(matrix_t matrix)
{
    int i, j;
    for(i = 0; i < matrix.n; i++) {
        for(j = 0; j < matrix.n; j++) {
            matrix.data[i][j] = (double)rand()/RAND_MAX*20 - 10;
        }
    }
}

/**
 *  Print the matrix.
 */
void matrix_print(matrix_t matrix)
{
    int i, j;
    for(i = 0; i < matrix.n; i++) {
        for(j = 0; j < matrix.n; j++) {
            printf("% .4f " , matrix.data[i][j]);
        }
        printf("\n");
    }
}

/**
 *  Compute the max norm of the matrix.
 */
double matrix_max_norm(matrix_t matrix)
{
    int i, j;
    double max = 0.0;
    max = matrix.data[0][0];
    for(i = 0; i < matrix.n; i++) {
        for(j = 0; j < matrix.n; j++) {
            if(matrix.data[i][j] > max) {
                max = matrix.data[i][j];
            }
        }
    }
    return max;
}


/**
* Matrix times Matrix
**/
matrix_t  matrix_mul(matrix_t A, matrix_t B,int n){
    matrix_t C ;
    int i = 0; int j = 0; int k = 0;
    matrix_allocate(&C,n); 
    matrix_fill_zeros(C);
    for (i = 0; i<n; i++){
        for(j =0; j<n; j++){
	    for (k = 0; k<n; k++)
            C.data[i][j] += A.data[i][k] * B.data[k][j];
        }
    }
    return C;
}

/**
 *  Fill the matrix with 0s.
 */
void matrix_fill_zeros(matrix_t matrix)
{
    int i, j;
    for(i = 0; i < matrix.n; i++) {
        for(j = 0; j < matrix.n; j++) {
            matrix.data[i][j] = 0;
        }
    }
}

