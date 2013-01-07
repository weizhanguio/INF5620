CC=cc
CFLAGS= -O2

FF=ftn
FFLAGS= -O2 

app: stream2_mpi.f
	$(CC) $(CGLAGS) -c mysecond.c
	$(FF) $(FFLAGS) -c stream2_mpi.f
	$(FF) $(FFLAGS)  stream2_mpi.o mysecond.o -o app

clean:
	rm -f app  *.o   *.e* 
