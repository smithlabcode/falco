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

EXTRA_DIST = README.md Configuration

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I $(top_srcdir)/src -DPROGRAM_PATH=\"$(abspath $(top_srcdir))\"

bin_PROGRAMS = falco

falco_CXXFLAGS = $(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)
falco_CPPFLAGS = -DPROGRAM_PATH=\"$(abspath $(top_srcdir))/Configuration\"
falco_SOURCES = \
	src/falco.cpp \
	src/HtmlMaker.cpp \
	src/FastqStats.cpp \
	src/StreamReader.cpp \
	src/FalcoConfig.cpp \
	src/OptionParser.cpp \
	src/smithlab_utils.cpp \
	src/HtmlMaker.hpp \
	src/FastqStats.hpp \
	src/StreamReader.hpp \
	src/FalcoConfig.hpp \
	src/OptionParser.hpp \
	src/smithlab_utils.hpp \
	src/aux.hpp
