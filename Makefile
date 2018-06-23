appname = friday
CC=clang++

LIBS=-lhts

SRC = $(shell find . -maxdepth 3 -name "*.cpp")
OBJ = $(patsubst %.cpp, %.o, $(SRC))

_DEPS = ./modules/handlers
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: clean $(appname)

$(appname): $(objects)
	$(CC) -std=c++11 -o $(appname) $^ $(LIBS) $(SRC)

.PHONY: clean

clean:
	rm -rf $(_DEPS)/*.o *.o $(appname)