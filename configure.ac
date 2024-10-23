dnl  falco: quality control for sequencing read files
dnl
dnl  Copyright (C) 2019-2024 Guilherme De Sena Brandine and
dnl                          Andrew D. Smith
dnl  Authors: Guilherme De Sena Brandine, Andrew Smith
dnl
dnl  This program is free software: you can redistribute it and/or
dnl  modify it under the terms of the GNU General Public License as
dnl  published by the Free Software Foundation, either version 3 of the
dnl  License, or (at your option) any later version.
dnl
dnl  This program is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl  General Public License for more details.

AC_INIT([falco], [1.2.4], [andrewds@usc.edu],
        [falco], [https://github.com/smithlabcode/falco])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17([noext], [mandatory])
AC_OPENMP([C++]) dnl make sure we have openmp for multi-core in falco

dnl check for the Zlib library
AC_CHECK_LIB([z], [zlibVersion], [],
  [AC_MSG_FAILURE([Zlib must be installed to build falco])], [])

dnl check for HTSLib if requested
hts_fail_msg="

Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"
AC_ARG_ENABLE([hts],
  [AS_HELP_STRING([--enable-hts], [Enable HTSLib @<:@yes@:>@])],
  [enable_hts=yes], [enable_hts=no])
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_LIB([hts], [hts_version], [],
  [AC_MSG_FAILURE([$hts_fail_msg])])]
)
AM_CONDITIONAL([ENABLE_HTS], [test "x$enable_hts" = "xyes"])

AC_CONFIG_FILES([Makefile])
dnl make the test data files available in the build tree
AC_CONFIG_LINKS([
test_build/md5sum.txt:test/md5sum.txt
test_build/test_data.tgz:test/test_data.tgz
])

AC_OUTPUT
