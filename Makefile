appname = friday
CC=clang++

LIBS=-lhts -lhdf5

SRC = $(shell find . -maxdepth 3 -name "*.cpp")
OBJ = $(patsubst %.cpp, %.o, $(SRC))

_DEPS = ./modules/handlers
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: clean $(appname)

$(appname): $(objects)
	$(CC) -std=c++11 -I./dep/hdf5-1.10.1/c++/src/ -I./dep/hdf5-1.10.1/build/hdf5/include/ -L./dep/hdf5-1.10.1/build/hdf5/lib/ -o $(appname) $^ $(LIBS) $(SRC)

.PHONY: clean

clean:
	rm -rf $(_DEPS)/*.o *.o $(appname)