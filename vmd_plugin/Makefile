CC ?= cc
LIB ?= toptraj.so
WARN := -Wall -Wextra -Werror -Wpedantic
CFLAGS += -fPIC -g -O2 -Itcl8.5
UNAME := $(shell uname -s)

LDFLAGS += -shared
ifeq ($(UNAME),Darwin)
LDFLAGS += -undefined dynamic_lookup
else
LDFLAGS += -Wl,-z,undefs
endif
LDLIBS += -lm

SRC := $(wildcard *.c)
HEADERS := $(wildcard *.h)

$(LIB): $(SRC) $(HEADERS)
	mkdir -p $(dir $@)
	$(CC) $(WARN) $(CFLAGS) $(SRC) -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f $(LIB)
