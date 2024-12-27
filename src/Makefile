CC = gcc
PKG_CFLAGS = -W -Wall -I$(R_HOME)/include -finline-functions -fPIC -std=gnu99 -O3
PKG_LIBS = -L$(R_HOME)/lib -lncurses -lpthread -lz -lm -lR

# OS := `uname`
# ifeq ($(OS),  Darwin)
# 	PKG_CFLAGS += -Wno-unused-function
# else
# 	PKG_LIBS += -lrt -ltinfo
# endif

SOURCES := $(wildcard *.c)
# SOURCES := bgzf.c, cfile.c, decompress.c, format.c, format2.c, format3.c, format6.c, format7.c, index.c, kstring.c, summary.c
OBJECTS := $(patsubst %.c, %.o, $(SOURCES))

knowYourCG.so: $(OBJECTS)
	$(CC) -shared -o $@ *.o $(PKG_LIBS)

# LHTSLIB_DIR = htslib
# LHTSLIB_INCLUDE = htslib/htslib
# LHTSLIB = $(LHTSLIB_DIR)/libhts.a
# $(LHTSLIB) :
# 	make -C $(LHTSLIB_DIR) libhts.a

%.o: %.c
	$(CC) -c $(PKG_CFLAGS) -I"." $< -o $@

# INCLUDE = include
# PROG = yame
# .PHONY: build
# build: exportcf $(PROG)

# exportcf:
# 	$(eval export CF_OPTIMIZE)

# .PHONY: debug
# debug: CF_OPTIMIZE := 0
# debug: CFLAGS += -g # -pg
# debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
# debug: build

#####################
##### libraries #####
#####################



# ###################
# ### subcommands ###
# ###################


# SOURCES := $(wildcard *.c)
# OBJECTS := $(patsubst %.c, %.o, $(SOURCES))

# LIBS=$(OBJECTS) $(LTHSLIB) # view.o chunk.o pack.o header.o bundle.o

# yame: $(LIBS)
# 	gcc $(CFLAGS) -o $@ *.o $(LTHSLIB) $(CLIB)


# ## clean just src
.PHONY: clean
clean :
	rm -f *.o *.so
## make -C $(LHTSLIB_DIR) clean
