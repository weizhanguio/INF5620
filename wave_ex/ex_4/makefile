LIBS=-lm
OBJS=wave1D_dn_ghost.o
CFLAGS=-O2
app:${OBJS}
	gcc -o app ${OBJS} ${LIBS} ${CFLAGS}
clean:
	rm -f app ${OBJS}

