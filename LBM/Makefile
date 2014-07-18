
CFLAGS1 = -g -O3
LIBS = -lm
CC	= g++
INCLUDES = globalParametersNoGraphics19.h

all: lb3d

clean:	
	rm *.o *.out

OBJS_1=   \
		LB_3D.o

lb3d:$(OBJS_1) 
	$(CC)  -o $@.out $(OBJS_1) $(CFLAGS1) $(LIBS)
%.o: %.c
	$(CC) -c $(FFLAGS1) $< -o $@ $(CFLAGS1) 

######### End of Makefile #########


