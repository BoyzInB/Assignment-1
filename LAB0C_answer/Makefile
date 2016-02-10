CC = gcc           # C compiler
LD = gcc           # C linker
CFLAGS =  -std=c99 -g -Wall  # Flags to be passed to the compiler
LDFLAGS =          # Flags for the linker
RM = /bin/rm -f    # Deletion program
OBJS = matrix.o matrix_main.o 
EXEC = matrix

all: matrix 

matrix: $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXEC) 

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c matrix.c

main.o: matrix.h matrix_main.c
	$(CC) $(CFLAGS) -c matrix_main.c
clean:
	$(RM) $(EXEC) $(OBJS)

