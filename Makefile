DIRS = -I/usr/local/include/ -L/usr/local/lib/ -Ithirdparty/emphf/ -Ithirdparty/streamvbyte/include/ -Lsdsl-inst/lib -Isdsl-inst/include
CFLAGS = $(DIRS) -O3 -march=native -g --std=c++11 -Wall -Wextra -Wno-implicit-fallthrough -Wno-ignored-qualifiers
LDFLAGS = $(DIRS) -lsdsl
CC = g++-9

CFILES = rmap.cpp lmer_reader.cpp # write_lmer.cpp mergetable.cpp index.cpp
HFILES = rmap.hpp lmer_reader.hpp write_lmer.hpp mergetable.hpp index.hpp
EXTRA = Makefile README.md COPYING
VERSION = 0.1

OBJS = $(CFILES:.cpp=.o)

all: build-index elmeri-dot candidates

build-index: $(OBJS) build-index.o
	$(CC) $(LDFLAGS) -o $@ $^ thirdparty/streamvbyte/streamvbyte_encode.o thirdparty/streamvbyte/streamvbyte_decode.o

candidates: $(OBJS) candidates.o
	$(CC) $(LDFLAGS) -o $@ $^ thirdparty/streamvbyte/streamvbyte_encode.o thirdparty/streamvbyte/streamvbyte_decode.o

elmeri-dot: $(OBJS) elmeri-dot.o
	$(CC) $(LDFLAGS) -o $@ $^ thirdparty/streamvbyte/streamvbyte_encode.o thirdparty/streamvbyte/streamvbyte_decode.o

%.o: %.cpp $(HFILES)
	$(CC) $(CFLAGS) -c $<

clean:
	rm $(OBJS)

dist:
	mkdir $(PROG)-$(VERSION)
	cp $(CFILES) $(PROG)-$(VERSION)/
	cp $(HFILES) $(PROG)-$(VERSION)/
	cp $(EXTRA) $(PROG)-$(VERSION)/
#	sed 's/VERSION/$(VERSION)/g' < README > $(PROG)-$(VERSION)/README
	tar zcvf $(PROG)-$(VERSION).tar.gz $(PROG)-$(VERSION)

dist-clean:
	rm -r $(PROG)-$(VERSION)/
