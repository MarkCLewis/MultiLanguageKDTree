# C Flags, -g, -Ofast
CFLAGS := -Wall -pedantic -std=c17 -g #-Ofast
CC := gcc
LDLIBS := -lm
INC :=

# run make -j 8, seems to be fastest

# ENDING is cpp

BUILD_DIR := build
SRCS := src/kdtree.c src/particle.c
#$(shell find src -iname '*.c' -type f)
HEADERS := $(shell find src -iname '*.h' -type f)
OBJS := $(SRCS:src/%.c=$(BUILD_DIR)/%.o)
OUT := kdtree-sim

MAIN_OBJ := $(BUILD_DIR)/kdtree-sim.o
TEST_OBJ := $(BUILD_DIR)/test.o

.PHONY: clean lines run


$(OUT): $(OBJS) $(MAIN_OBJ)
	$(CC) -o $@ $^ $(LDLIBS)

# run-test: test
# 	./$<

test: $(OBJS) $(TEST_OBJ)
	$(CC) -o $@ $^ $(LDLIBS)

$(BUILD_DIR)/%.o: src/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) -c $(INC) -o $@ $< $(CFLAGS)

clean:
	rm -rf $(BUILD_DIR) $(OUT)

lines: $(SRCS) $(HEADERS)
	wc -l $^ | grep total
