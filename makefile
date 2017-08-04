CC=g++
CFLAGS=-I. -g -O0
CLIBS=-lm
DEPS = XGridUtil.h
OBJ = XGridUtil.o testxgrid.o

%.o: %.C $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

testxgrid: $(OBJ)
	$(CC) -o $@ $^ $(CLIBS)

clean:
	rm *.o testxgrid
