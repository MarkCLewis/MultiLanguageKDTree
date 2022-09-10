# C Flags, -g, -Ofast
CFLAGS := -Wall -pedantic -std=c17 -Ofast
CC := gcc
LDLIBS := -lm
INC :=

# run make -j 8, seems to be fastest

# ENDING is cpp

BUILD_DIR := build
SRCS := $(shell find src -iname '*.c' -type f)
HEADERS := $(shell find src -iname '*.h' -type f)
OBJS := $(SRCS:src/%.c=$(BUILD_DIR)/%.o)
OUT := kdtree-sim

.PHONY: clean lines run


$(OUT): $(OBJS)
	$(CC) -o $@ $^ $(LDLIBS)

$(BUILD_DIR)/%.o: src/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) -c $(INC) -o $@ $< $(CFLAGS)

clean:
	rm -rf $(BUILD_DIR) $(OUT)

lines: $(SRCS) $(HEADERS)
	wc -l $^ | grep total
