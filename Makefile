INCDIRS = -I/usr/local/include/ -Ithirdparty/emphf/ -Ithirdparty/streamvbyte/include/ -Isdsl-inst/include
LINKDIRS = -L/usr/local/lib/ -Lsdsl-inst/lib
OPT = -O3 -march=native -msse4.2
CFLAGS = $(INCDIRS) -g --std=c++11 -Wall -Wextra -Wno-implicit-fallthrough -Wno-ignored-qualifiers $(OPT)
LDFLAGS = $(LINKDIRS) -lsdsl thirdparty/streamvbyte/streamvbyte_encode.o thirdparty/streamvbyte/streamvbyte_decode.o
CXX = g++

CFILES = rmap.cpp lmer_reader.cpp
HFILES = rmap.hpp lmer_reader.hpp write_lmer.hpp mergetable.hpp index.hpp
EXTRA = Makefile README.md COPYING
VERSION = 1.0

OBJS = $(CFILES:.cpp=.o)

all: selkie-index selkie-dot selkie-candidates

selkie-index: $(OBJS) build-index.o
	$(CXX) $(LDFLAGS) -o $@ $^

selkie-candidates: $(OBJS) candidates.o
	$(CXX) $(LDFLAGS) -o $@ $^

selkie-dot: $(OBJS) dot.o
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp $(HFILES)
	$(CXX) $(CFLAGS) -c $<

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
