FLAGS := -Wall -pedantic -std=c++17 -Ofast -fomit-frame-pointer
SRC_FILES := $(wildcard *.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))
CC := clang++

kdtree-sim: $(OBJ_FILES)
	$(CC) $(FLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(FLAGS) -c -o $@ $<

clean:
	rm *.o
	rm kdtree-sim

all:
	make kdtree-sim
