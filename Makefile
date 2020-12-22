# Copyright (C) 2019-2020 Guilherme de Sena Brandine
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
#

SRC_ROOT=$(shell pwd)
all:
	@make -C src SRC_ROOT=$(SRC_ROOT) all

install:
	@make -C src SRC_ROOT=$(SRC_ROOT) install

clean:
	@make -C src clean
.PHONY: clean
