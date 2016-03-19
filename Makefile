CC = g++

CLASSDIR = /home/cristinel/noc/remnoc
INCDIRS = $(CLASSDIR)/include

LIB_DIR = -L/usr/lib/X11
LIB = -lX11 -lm
X11_INCLUDE = -I/usr/X11R6/include

WARN_FLAGS = -Wall -Wpointer-arith -Wcast-qual -Wstrict-prototypes -O -D__USE_FIXED_PROTOTYPES__ -ansi -pedantic -Wmissing-prototypes -Wshadow -Wcast-align -D_POSIX_SOURCE
DEBUG_FLAGS = -g
OPT_FLAGS = -O3

FLAGS = $(OPT_FLAGS)
FLAGS += $(addprefix -I, $(INCDIRS))

EXE = remnoc

OBJ = remnoc_hungarian.o remnoc_utils.o remnoc.o remnoc_annealer.o remnoc_main.o remnoc_gui.o 

SRC = remnoc_hungarian.cpp remnoc_utils.cpp remnoc.cpp remnoc_annealer.cpp remnoc_main.cpp remnoc_gui.cpp

H = include/remnoc_hungarian.h include/remnoc_utils.h include/remnoc.h include/remnoc_annealer.h include/remnoc_gui.h 

$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE) $(LIB_DIR) $(LIB)

remnoc_hungarian.o: remnoc_hungarian.cpp $(H)
	$(CC) -c $(FLAGS) remnoc_hungarian.cpp

remnoc_utils.o: remnoc_utils.cpp $(H)
	$(CC) -c $(FLAGS) remnoc_utils.cpp

remnoc.o: remnoc.cpp $(H)
	$(CC) -c $(FLAGS) remnoc.cpp

remnoc_annealer.o: remnoc_annealer.cpp $(H)
	$(CC) -c $(FLAGS) remnoc_annealer.cpp

remnoc_main.o: remnoc_main.cpp $(H)
	$(CC) -c $(FLAGS) remnoc_main.cpp

remnoc_gui.o: remnoc_gui.cpp $(H)
	$(CC) -c $(FLAGS) $(X11_INCLUDE) remnoc_gui.cpp
