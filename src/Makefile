# falco: quality control for sequencing read files
#
# Copyright (C) 2019 Guilherme De Sena Brandine and
#                    Andrew D. Smith
# Authors: Guilherme De Sena Brandine, Andrew Smith
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

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
curr_makefile := $(abspath $(lastword $(MAKEFILE_LIST)))
current_directory := $(dir $(curr_makefile))/..
CXXFLAGS += -DPROGRAM_PATH=\"$(current_directory)\"

all: $(PROGRAMS)

install: $(PROGRAMS)
	@mkdir -p bin
	@install -m 755 $(PROGRAMS) bin
	rm $(PROGRAMS)

$(PROGRAMS): HtmlMaker.o FastqStats.o StreamReader.o FalcoConfig.o \
	OptionParser.o smithlab_utils.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	@-rm -f $(PROGRAMS) *.o
.PHONY: clean
