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

EXTRA_DIST = README.md LICENSE Configuration example.fq test

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I $(top_srcdir)/src

### ADS: use AM_CXXFLAGS if these are needed
## AM_CXXFLAGS = -Wall -Wextra -Wpedantic -Wno-unknown-attributes

TESTS = test/falco.test
TEST_EXTENSIONS = .test

bin_PROGRAMS = falco

falco_CXXFLAGS = $(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)
falco_CPPFLAGS = -DPROGRAM_PATH=\"$(abspath $(top_srcdir))\"
if ENABLE_HTS
falco_CPPFLAGS += -DUSE_HTS
endif

falco_SOURCES = \
	src/falco.cpp \
	src/FastqStats.cpp \
	src/HtmlMaker.cpp \
	src/Module.cpp \
	src/StreamReader.cpp \
	src/FalcoConfig.cpp \
	src/OptionParser.cpp \
	src/smithlab_utils.cpp \
	src/Module.hpp \
	src/FastqStats.hpp \
	src/HtmlMaker.hpp \
	src/StreamReader.hpp \
	src/FalcoConfig.hpp \
	src/OptionParser.hpp \
	src/smithlab_utils.hpp \
	src/aux.hpp

CLEANFILES = tests_build/test_data.tgz
