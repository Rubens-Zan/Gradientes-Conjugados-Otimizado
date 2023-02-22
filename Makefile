# Autor: Rubens Laszlo
# Data: 12/2022
# GRR 20206147 

CC = gcc
EXEC = cgSolver
CFLAG = -Wall -std=c99 -lm -O3 -mavx -march=native -DLIKWID_PERFMON -I/home/soft/likwid/include  
MODULOS = sislin \
	utils \
	resolvedorGradConjug 
	
OBJETOS = main.o $(addsuffix .o,$(MODULOS))

.PHONY: all clean debug


all: CGSOLVER

debug: CFLAGS += -DDEBUG
debug: CGSOLVER

CGSOLVER: $(OBJETOS)
	$(CC) -o $(EXEC) $(OBJETOS) $(CFLAG)

clean:
	$(RM) $(OBJETOS) $(EXEC)