CFLAGS =  -O3 -fPIE
# -ggdb
LINKLIBS= -lm


OBJ0 = RngStream.o gaussian-r-n.o array_alloc.o vectors2d.o useful-tools.o sim-gauss-2.o

G-2:    $(OBJ0)
	gcc -o $@ $(OBJ0) $(LINKLIBS)
	
clean:
	rm -f $(OBJ0) G-2

