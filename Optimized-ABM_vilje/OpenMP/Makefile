CC=icc


app:  temporalblocking_omp.c
	$(CC)  temporalblocking_omp.c  -o app -vec-report2 -openmp -O2

clean:
	rm -rf app  *.o  *.*~
