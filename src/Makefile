CC = gcc
PKG_CFLAGS = -W -Wall -I$(R_HOME)/include -finline-functions -fPIC -std=gnu99 -O3
PKG_LIBS = -L$(R_HOME)/lib -lncurses -lpthread -lz -lm -lR

SOURCES := $(wildcard *.c)
OBJECTS := $(patsubst %.c, %.o, $(SOURCES))
# SOURCES := bgzf.c, cfile.c, decompress.c, format.c, format2.c, format3.c, format6.c, format7.c, index.c, kstring.c, summary.c

knowYourCG.so: $(OBJECTS)
	$(CC) -shared -o $@ *.o $(PKG_LIBS)

%.o: %.c
	$(CC) -c $(PKG_CFLAGS) -I"." $< -o $@

# ## clean just src
.PHONY: clean
clean :
	rm -f *.o *.so
