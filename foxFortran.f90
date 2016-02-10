PROGRAM FOX
IMPLICIT NONE
include "mpif.h"

INTEGER :: procs, rank, error
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

INTEGER :: g_order, g_side, my_g_row, my_g_column, my_g_rank, comm, row_comm, col_comm, block_mpi_t
INTEGER :: rows
REAL, DIMENSION(:,:), ALLOCATABLE :: matrixA, matrixB, matrixC, localA, tempA, localB, localC

CALL MPI_Init ( error )
CALL MPI_Comm_size ( MPI_COMM_WORLD, procs, error )
CALL MPI_Comm_rank ( MPI_COMM_WORLD, rank, error )

CALL Read_Matrix(rank,error,rows,matrixA,matrixB,matrixC)
CALL Setup_grid(g_order,g_side,error,comm,my_g_rank,my_g_row,my_g_column,row_comm,col_comm)
CALL Distribute_matrixes(g_side,rows,block_mpi_t,error,matrixA,matrixB,localA,localB,localC,my_g_row,my_g_column)
CALL Perform_Fox_Algorithm(my_g_row,my_g_column,g_order,localA,tempA,localB,localC,row_comm,col_comm,status,error)
CALL Print_Last_Result(rank, my_g_rank, g_side, comm, procs, error, status, localC, matrixC, tempA)


CALL MPI_Finalize (error)

CONTAINS

SUBROUTINE Print_Last_Result(rank, my_g_rank, g_side, comm, procs, error, status, localC, matrixC, tempA)
INTEGER, INTENT(IN) :: rank, my_g_rank, g_side, comm, procs
INTEGER, INTENT(OUT) :: error
INTEGER, DIMENSION(:), INTENT(OUT) :: status
REAL, DIMENSION(:,:), INTENT(IN) :: localC
REAL, DIMENSION(:,:), INTENT(OUT) :: matrixC, tempA

INTEGER :: grow,gcol,i
INTEGER,DIMENSION(2) :: coordinates

IF (rank .EQ. 0 .AND. my_g_rank .NE. 0) PRINT*, "Houston, we have a problem with I/O"

CALL MPI_SEND(localC,g_side*g_side,MPI_REAL,0,0,comm,error)

IF (rank .EQ. 0) THEN
matrixC=0
DO i=1,procs
CALL MPI_Recv(tempA,g_side*g_side,MPI_REAL,MPI_ANY_SOURCE,0,comm,status,error)
CALL MPI_Cart_coords(comm,status(MPI_Source),2,coordinates,error)
grow = coordinates(1)
gcol = coordinates(2)
matrixC(grow*g_side+1:(grow+1)*g_side,gcol*g_side+1:(gcol+1)*g_side) = tempA
END DO
PRINT*,"============ Matrix multiplication in parallel ==================="
CALL Print_Matrix(matrixC)
END IF
END SUBROUTINE Print_Last_Result


SUBROUTINE Print_Matrix(M)
REAL, DIMENSION(:,:) :: M
INTEGER :: side,row,column

side = SIZE(M,1)

DO row=1,side
DO column=1,side-1
WRITE(*,'(F10.3)',ADVANCE='NO'), M(row,column)
END DO
WRITE(*,'(F10.3)'), M(row,side)
END DO
END SUBROUTINE Print_Matrix


SUBROUTINE Perform_Fox_Algorithm(my_g_row,my_g_column,g_order,localA,tempA,localB,localC,row_comm,col_comm,status,error)
INTEGER, INTENT(IN) :: my_g_row,my_g_column,g_order,row_comm,col_comm
INTEGER, INTENT(OUT) :: error
INTEGER, DIMENSION(:), INTENT(OUT) :: status
REAL, DIMENSION(:,:), INTENT(IN) :: localA
REAL, DIMENSION(:,:), INTENT(INOUT) :: tempA,localB,localC

INTEGER :: source,destination,stage,bcast_root

localC = 0
source = MOD(my_g_row + 1,g_order)
destination = MOD(my_g_row + g_order - 1,g_order)

DO stage=0,g_order-1
bcast_root = MOD(my_g_row + stage,g_order)

IF (my_g_column .EQ. bcast_root) THEN
CALL MPI_BCAST(localA,g_side*g_side,MPI_REAL,bcast_root,row_comm,error)
CALL Matrix_Multiply(localA,localB,localC)
ELSE
CALL MPI_BCAST(tempA,g_side*g_side,MPI_REAL,bcast_root,row_comm,error)
CALL Matrix_Multiply(tempA,localB,localC)
END IF
CALL MPI_Sendrecv_replace(localB,g_side*g_side,MPI_REAL,destination,0,source,0,col_comm,status,error)
END DO
END SUBROUTINE Perform_Fox_Algorithm


SUBROUTINE Matrix_Multiply(A,B,C)
REAL, DIMENSION(:,:), INTENT(IN) :: A,B
REAL, DIMENSION(:,:), INTENT(OUT) :: C
INTEGER :: side,row,column

side = SIZE(A,1)

DO row=1,side
DO column=1,side
C(row,column) = C(row,column) + SUM(A(row,:)*B(:,column))
END DO
END DO
END SUBROUTINE Matrix_Multiply


SUBROUTINE Distribute_matrixes(g_side,rows,block_mpi_t,error,matrixA,matrixB,localA,localB,localC,my_g_row,my_g_column)
INTEGER, INTENT(IN) :: g_side,rows,my_g_row,my_g_column
INTEGER, INTENT(OUT) :: block_mpi_t,error
REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: localA,localB,localC
REAL, DIMENSION(:,:), INTENT(INOUT) :: matrixA,matrixB

CALL MPI_TYPE_VECTOR(g_side,g_side,rows,MPI_REAL,block_mpi_t,error)
CALL MPI_TYPE_COMMIT(block_mpi_t, error)

CALL MPI_BCAST(matrixA,rows*rows,MPI_REAL,0,MPI_COMM_WORLD,error)
CALL MPI_BCAST(matrixB,rows*rows,MPI_REAL,0,MPI_COMM_WORLD,error)

ALLOCATE(localA(g_side,g_side),tempA(g_side,g_side),localB(g_side,g_side),localC(g_side,g_side))
localA = matrixA(my_g_row*g_side+1:(my_g_row+1)*g_side,my_g_column*g_side+1:(my_g_column+1)*g_side)
localB = matrixB(my_g_row*g_side+1:(my_g_row+1)*g_side,my_g_column*g_side+1:(my_g_column+1)*g_side)
END SUBROUTINE Distribute_matrixes


SUBROUTINE Setup_grid(g_order,g_side,error,comm,my_g_rank,my_g_row,my_g_column,row_comm,col_comm)
INTEGER, INTENT(OUT) :: g_order, g_side,error,comm,my_g_rank,my_g_row,my_g_column,row_comm,col_comm

INTEGER,DIMENSION(2) :: dimensions,coordinates
LOGICAL,DIMENSION(2) :: wrap_around,free_coords

g_order = SQRT(REAL(procs))
g_side = rows / g_order

dimensions = g_order
wrap_around = .TRUE.
CALL MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_around,.TRUE.,comm,error)
CALL MPI_Comm_rank ( comm, my_g_rank, error )
CALL MPI_Cart_coords(comm,my_g_rank,2,coordinates,error)
my_g_row = coordinates(1)
my_g_column = coordinates(2)

free_coords(1) = .FALSE.
free_coords(2) = .TRUE.
CALL MPI_Cart_sub(comm,free_coords,row_comm,error)

free_coords(1) = .TRUE.
free_coords(2) = .FALSE.
CALL MPI_Cart_sub(comm,free_coords,col_comm,error)
END SUBROUTINE Setup_grid

SUBROUTINE Read_Matrix(rank,error,rows,matrixA,matrixB,matrixC)
INTEGER, INTENT(IN) :: rank
INTEGER, INTENT(OUT) :: rows,error
INTEGER :: i
REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: matrixA,matrixB,matrixC

IF (rank .EQ. 0) THEN
READ*, rows
CALL MPI_BCAST(rows,1,MPI_REAL,0,MPI_COMM_WORLD,error)
ALLOCATE(matrixA(rows,rows),matrixB(rows,rows),matrixC(rows,rows))
DO i=1,rows
READ*, matrixA(i,:)
END DO
DO i=1,rows
READ*, matrixB(i,:)
END DO

!! Calculate the matrix multiplication for comparison
matrixC = 0
CALL Matrix_Multiply(matrixA,matrixB,matrixC)

PRINT*, "CALCULATED IN PROCESS 0"
PRINT*, "=======MATRIX A=========="
CALL Print_Matrix(matrixA)
PRINT*, "=======MATRIX B=========="
CALL Print_Matrix(matrixB)
PRINT*, "=======MATRIX C=========="
CALL Print_Matrix(matrixC)
ELSE
CALL MPI_BCAST(rows,1,MPI_REAL,0,MPI_COMM_WORLD,error)
ALLOCATE(matrixA(rows,rows),matrixB(rows,rows))
END IF
END SUBROUTINE Read_Matrix

END PROGRAM FOX