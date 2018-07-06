appname = friday
CC=g++-8

LIBS=-lhts -lhdf5

SRC = $(shell find . -maxdepth 3 -name "*.cpp" -not -path "./dep/*")
OBJ = $(patsubst %.cpp, %.o, $(SRC))

_DEPS = ./modules/handlers
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: clean $(appname)

$(appname): $(objects)
	$(CC) -std=c++11 -I./dep/htslib-1.8/include/ -I./dep/cpptqdm/ -I./dep/hdf5-1.10.2/build/hdf5/include/ -L./dep/htslib-1.8/lib/ -L./dep/hdf5-1.10.2/build/hdf5/lib/ -o $(appname) $^ $(LIBS) $(SRC) -fopenmp

.PHONY: clean

clean:
	rm -rf $(_DEPS)/*.o *.o ./outputs/ $(appname)
	mkdir ./outputs/