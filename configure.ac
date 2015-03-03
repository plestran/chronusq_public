#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT([Chrous Quantum Chemistry], [1.0], [xsli@uw.edu])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE
AC_CONFIG_SRCDIR([include/aointegrals.h])
AC_CONFIG_HEADERS([config.h])

: ${LIBS+="-llapack -lblas"}
#: ${LIBS="-lm -lc -pgf77libs"}
#: ${CFLAGS="-I../../include -w "}
#: ${CXXFLAGS="-I../../include -w "}
#: ${F77FLAGS=""}

#AM_CONDITIONAL([PGI],[test "x$F77"=xpgf77])
#AM_COND_IF([PGI],[CFLAGS="-w -mp=nonuma  -tp nehalem-64 -Mreentrant -Minfo -Mneginfo -time -fast -O2 -Munroll -Mvect=assoc,recog,cachesize:20971520"],[])
#AM_COND_IF([PGI],[CXXFLAGS="-w -mp=nonuma  -tp nehalem-64 -Mreentrant -Minfo -Mneginfo -time -fast -O2 -Munroll -Mvect=assoc,recog,cachesize:20971520"],[])
#AM_COND_IF([PGI],[LIBS+="-lm -lc -pgf77libs"],[])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([blas], [dgemm_])
AC_CHECK_LIB([lapack], [dsyev_])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/aointegrals/Makefile
                 src/atlas/Makefile
                 src/basisset/Makefile
                 src/cerr/Makefile
                 src/chronQ/Makefile
                 src/classtools/Makefile
                 src/controls/Makefile
                 src/fileio/Makefile
                 src/matrix/Makefile
                 src/mattest/Makefile
                 src/memory/Makefile
                 src/mointegrals/Makefile
                 src/molecule/Makefile
                 src/sdresponse/Makefile
                 src/singleslater/Makefile
                 src/tools/Makefile
                 src/worker/Makefile])
AC_OUTPUT
