CC=cc
CFLAGS=  -O2

FF=ftn
FFLAGS= -O2

app: stream2.f
	$(CC) $(CFLAGS)  -c mysecond.c
	$(FF) $(FFLAGS) -c stream2.f
	$(FF) $(FFLAGS)  stream2.o mysecond.o -o app

clean:
	rm -f app  *.o    *.e*
