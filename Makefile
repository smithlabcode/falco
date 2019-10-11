PROGRAMS = falco
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++11

# optional zlib and hts
LIBS =
ifndef NO_ZLIB
LIBS += -lz
CXXFLAGS += -DUSE_ZLIB
endif

ifndef NO_HTS
LIBS += -lhts
CXXFLAGS += -DUSE_HTS
endif

# get path to binary at compile time to look for files
CXXFLAGS += "-DPROGRAM_PATH=\"`pwd`\""

all: $(PROGRAMS)

install: $(PROGRAMS)
	@mkdir -p bin
	@install -m 755 $(PROGRAMS) bin
	rm $(PROGRAMS)

$(PROGRAMS): $(addprefix src/, HtmlMaker.o \
	FastqStats.o StreamReader.o Config.o OptionParser.o smithlab_utils.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -rf bin/$(PROGRAMS) src/*.o

.PHONY: clean
