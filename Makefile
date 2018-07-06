appname = friday.exe
CC=g++-8
INC=-I./dep/htslib-1.8/include/ -I./dep/cpptqdm/ -I./dep/hdf5-1.10.2/build/hdf5/include/
LIBS=-lhts -lhdf5 -fopenmp -L./dep/htslib-1.8/lib/ -L./dep/hdf5-1.10.2/build/hdf5/lib/

SRC = $(shell find . -maxdepth 3 -name "*.cpp" -not -path "./dep/*")
OBJ = $(patsubst %.cpp, %.o, $(SRC))

_DEPS = ./modules/handlers
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: clean $(appname)

$(appname): $(objects)
	$(CC) -std=c++11 $(INC) -o $(appname) $^ $(SRC) $(LIBS)

.PHONY: clean
.PHONY: makedeps
.PHONY: cleandeps

makedeps: cleandeps
	cd ./dep/ && \
	tar -xvzf hdf5-1.10.2.tar.gz && \
	cd hdf5-1.10.2 && \
	mkdir build && \
	cd ./build/ && \
	../configure && \
	make && \
	make install

	rm -rf ./dep/htslib-1.8
	cd ./dep/ && \
    tar -xvjf htslib-1.8.tar.bz2 && \
    cd htslib-1.8 && \
    ./configure --prefix=$(shell pwd)/dep/htslib-1.8/ && \
    make && \
    make install

cleandeps:
	rm -rf ./dep/htslib-1.8
	rm -rf ./dep/hdf5-1.10.2

clean:
	rm -rf $(_DEPS)/*.o *.o ./outputs/ $(appname)
	mkdir ./outputs/