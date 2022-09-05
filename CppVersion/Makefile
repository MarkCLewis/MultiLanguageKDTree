FLAGS := -Wall -pedantic -std=c++17 -Ofast
SRC_FILES := $(wildcard *.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(SRC_FILES))

kdtree-sim: $(OBJ_FILES)
	g++ $(FLAGS) -o $@ $^

%.o: %.cpp
	g++ $(FLAGS) -c -o $@ $<

clean:
	rm *.o
	rm kdtree-sim

all:
	make kdtree-sim
