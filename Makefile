INCLUDE_DIR = ./include
CFLAGS = -std=c11 -I$(INCLUDE_DIR) #-Wall
VPATH = ./src

CC = gcc

TARGETS = 3BP

DEPENDENCIES = main.c main_functions.c auxiliar_functions.c auxiliar_functions_gsl.c auxiliar_functions_vmo.c dynamical_system.c

build: $(TARGETS)

$(TARGETS): $(DEPENDENCIES) -lgsl -lgslcblas -lm  
			$(CC) $(CFLAGS) -o $(TARGETS) $^

.PHONY: clean
clean:
	-rm -f $(TARGETS)
