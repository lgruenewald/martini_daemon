CC = gcc
TCL_VERSION := 8.6
CFLAGS := -shared $(shell pkg-config --cflags "tcl$(TCL_VERSION)") -fPIC -g -O1 -Wall -Wextra
LDFLAGS := $(shell pkg-config --libs "tcl$(TCL_VERSION)")
LIB := toptraj.so

SRC := $(wildcard *.c)
HEADERS := $(wildcard *.h)

.PHONY: all
all: lib

.PHONY: lib
lib: $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(LIB) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(LIB)
