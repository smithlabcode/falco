# Copyright (C) 2018-2019 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This file is part of ABISMAL.
#
# ABISMAL is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ABISMAL is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.

PROGS = falco

CXX = g++
CXXFLAGS = -Wall -std=c++11
CPPFLAGS = -DPROGRAM_PATH=\"$(SRC_ROOT)\"
LDLIBS = -lz
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef HAVE_HTSLIB
CPPFLAGS += -DUSE_HTS
LDLIBS += -lhts
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)
install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS)

$(PROGS): FalcoConfig.o FastqStats.o HtmlMaker.o Module.o OptionParser.o \
	        smithlab_utils.o StreamReader.o

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS) $(LDLIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~
.PHONY: clean

