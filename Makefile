DIRS = -I/usr/local/include/ -L/usr/local/lib/ -Ithirdparty/emphf/ -Ithirdparty/streamvbyte/include/ -Lsdsl-inst/lib -Isdsl-inst/include
CFLAGS = $(DIRS) -O3 -march=native -g --std=c++11 -Wall -Wextra -Wno-implicit-fallthrough -Wno-ignored-qualifiers
LDFLAGS = $(DIRS) -lsdsl thirdparty/streamvbyte/streamvbyte_encode.o thirdparty/streamvbyte/streamvbyte_decode.o
CC = g++-9

CFILES = rmap.cpp lmer_reader.cpp
HFILES = rmap.hpp lmer_reader.hpp write_lmer.hpp mergetable.hpp index.hpp
EXTRA = Makefile README.md COPYING
VERSION = 1.0

OBJS = $(CFILES:.cpp=.o)

all: selkie-index selkie-dot selkie-candidates

selkie-index: $(OBJS) build-index.o
	$(CC) $(LDFLAGS) -o $@ $^

selkie-candidates: $(OBJS) candidates.o
	$(CC) $(LDFLAGS) -o $@ $^

selkie-dot: $(OBJS) dot.o
	$(CC) $(LDFLAGS) -o $@ $^

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
