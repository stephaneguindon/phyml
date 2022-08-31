#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([PhyML],esyscmd([sh -c "echo "3.3." | tr -d '\n' ; echo "20220408" | tr -d '\n'"]),[guindon@lirmm.fr])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/main.c],[doc/phyml-manual.tex])
AC_CONFIG_HEADERS([config.h])

AC_DEFINE([UNIX],[1],[Unix tag on])
AC_DEFINE([DEBUG],[1],[Debug tag on])

AC_CANONICAL_HOST

AC_PROG_CC([gcc clang mpicc])

AS_IF([test "$CC" = "gcc"],[DEFAULT_VECTOR_FLAG="-march=native"])
AS_IF([test "$CC" = "clang"],[DEFAULT_VECTOR_FLAG="-mcpu=native"])
AS_IF([test "$CC" = "mpicc"],[DEFAULT_VECTOR_FLAG="-march=native"])

AC_MSG_NOTICE(["$CC"])


# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Check for presence of LaTeX and pdfLaTeX
AC_CHECK_PROGS([PDFLATEX], pdflatex, [/bin/false])
AC_CHECK_PROGS([LATEX], latex, [/bin/false])
if test "x$PDFLATEX" = x/bin/false || test "x$LATEX" = x/bin/false; then
 enable_docs=no
AC_MSG_WARN([Documentation will not be built!])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])


dnl DEFAULT_VECTOR_FLAG=""
dnl DEFAULT_VECTOR_FLAG="-mavx -mfma"
dnl DEFAULT_VECTOR_FLAG="-mavx2"
dnl DEFAULT_VECTOR_FLAG=-msse3
dnl DEFAULT_VECTOR_FLAG="-march=native"
dnl DEFAULT_VECTOR_FLAG="-mtune=native"
dnl DEFAULT_VECTOR_FLAG="-mcpu=native"

AC_ARG_ENABLE([native],
              [AS_HELP_STRING([--disable-native],
                              [Remove AVX or SSE compilation options.])])

AS_IF([test "x$enable_native" = "xno"],
      [VECTOR_FLAG="-DDISABLE_NATIVE"],
      [VECTOR_FLAG="${DEFAULT_VECTOR_FLAG}"])


LDFLAGS="${ARGC_flag}" 

dnl Uncomment two rows below when compiling rwarp (need to have libtool installed)
dnl LT_INIT
dnl AC_PROG_LIBTOOL


dnl DEFAULT_C_FLAG="-fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls -O2 -fno-omit-frame-pointer -g ${VECTOR_FLAG}"
DEFAULT_C_FLAG="-std=c99 -O3 -fomit-frame-pointer -funroll-loops -Wall -Winline -finline ${VECTOR_FLAG}"
dnl DEFAULT_C_FLAG="-O2 ${VECTOR_FLAG}"
dnl DEFAULT_C_FLAG="${VECTOR_FLAG}"
dnl DEFAULT_C_FLAG="-pthread  -O3 -g -ffunction-sections -fdata-sections"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Remove optimization options and add debug informations.])])
AS_IF([test "x$enable_debug" = "xyes"],
      [CFLAGS="-ansi ${VECTOR_FLAG} -pedantic -Wall -std=c99 -O0 -g -Wuninitialized -Wunused -Wunreachable-code"],
      [CFLAGS="${DEFAULT_C_FLAG}"])

AC_ARG_ENABLE([safemode],
        AS_HELP_STRING([--enable-safemode],[Enable safe mode.]))
AS_IF([test "x$enable_safemode" = "xyes"],
      [CFLAGS="-Wall -DSAFEMODE -O3 ${VECTOR_FLAG}"],
      [AS_IF([test "x$enable_debug" = "xno"],[CFLAGS="${DEFAULT_C_FLAG}"])])

AC_ARG_ENABLE([gprof],
              [AS_HELP_STRING([--enable-gprof],
                              [Remove optimization options and add profiling informations.])])
AS_IF([test "x$enable_gprof" = "xyes"],
      [CFLAGS="-g -pg -fno-inline -fprofile-arcs -ftest-coverage ${VECTOR_FLAG}"],
      [AS_IF([test "x$enable_debug" = "xno"],[CFLAGS="${DEFAULT_C_FLAG}"])])

AC_ARG_ENABLE([phyml-mpi],
              [AS_HELP_STRING([--enable-phyml-mpi],
                              [Compile with mpicc instead of gcc.])])
AS_IF([test "x$enable_phyml_mpi" = "xyes"],[CC="mpicc"])
AS_IF([test "x$enable_phyml_mpi" = "xyes"],AC_DEFINE([MPI],[1],[MPI tag on]))
AM_CONDITIONAL([WANT_MPI], [test "$enable_phyml_mpi" = yes])


AC_ARG_ENABLE([win],
              [AS_HELP_STRING([--enable-win],
                              [Compile with mingw instead of gcc.])])
AS_IF([test "x$enable_win" = "xyes"],[CC="i686-w64-mingw32-gcc"])
AS_IF([test "x$enable_win" = "xyes"],AC_DEFINE([WIN32],[1],[WIN32 tag on]))
AM_CONDITIONAL([WANT_WIN], [test "$enable_win" = yes])

AC_ARG_ENABLE([phytime],[AS_HELP_STRING([--enable-phytime],[Compile PhyTime])],[phytime=yes],[phytime=no])
AM_CONDITIONAL([WANT_PHYTIME], [test "$phytime" = yes])
if test "$phytime" = yes; then
   AC_DEFINE([PHYTIME],[1],[PHYTIME tag on])
fi

AC_ARG_ENABLE([phyml],[AS_HELP_STRING([--enable-phyml],[Compile PhyML])],[phyml=yes],[phyml=no])
AM_CONDITIONAL([WANT_PHYML], [test "$phyml" = yes])
if test "$phyml" = yes; then
   AC_DEFINE([PHYML],[1],[PHYML tag on])
fi


AC_ARG_ENABLE([rwrap],[AS_HELP_STRING([--enable-rwrap],[Compile Rwrap])],[rwrap=yes],[rwrap=no])
AM_CONDITIONAL([WANT_RWRAP], [test "$rwrap" = yes])
if test "$rwrap" = yes; then
   AC_DEFINE([RWRAP],[1],[RWRAP tag on])
fi

AC_ARG_ENABLE([rf],[AS_HELP_STRING([--enable-rf],[Compile RF])],[rf=yes],[rf=no])
AM_CONDITIONAL([WANT_RF], [test "$rf" = yes])
if test "$rf" = yes; then
   AC_DEFINE([RF],[1],[RF tag on])
fi

AC_ARG_ENABLE([test],[AS_HELP_STRING([--enable-test],[Compile test])],[test=yes],[test=no])
AM_CONDITIONAL([WANT_TEST], [test "$test" = yes])
if test "$test" = yes; then
   AC_DEFINE([TEST],[1],[TEST tag on])
fi

AC_ARG_ENABLE([evolve],[AS_HELP_STRING([--enable-evolve],[Compile evolve])],[evolve=yes],[evolve=no])
AM_CONDITIONAL([WANT_EVOLVE], [test "$evolve" = yes])
if test "$evolve" = yes; then
   AC_DEFINE([EVOLVE],[1],[EVOLVE tag on])
fi

AC_ARG_ENABLE([phyrex],[AS_HELP_STRING([--enable-phyrex],[Compile phyrex])],[phyrex=yes],[phyrex=no])
AM_CONDITIONAL([WANT_PHYREX], [test "$phyrex" = yes])
if test "$phyrex" = yes; then
   AC_DEFINE([PHYREX],[1],[PHYREX tag on])
fi

AC_ARG_ENABLE([phyrexsim],[AS_HELP_STRING([--enable-phyrexsim],[Compile phyrexsim])],[phyrexsim=yes],[phyrexsim=no])
AM_CONDITIONAL([WANT_PHYREXSIM], [test "$phyrexsim" = yes])
if test "$phyrexsim" = yes; then
   AC_DEFINE([PHYREXSIM],[1],[PHYREXSIM tag on])
fi

dnl AS_IF([test "x$enable_phyrex" = "xyes"],[PKG_CHECK_MODULES([GTK], [gtk+-3.0])])
dnl AS_IF([test "x$enable_phyrex" = "xyes"],[PKG_CHECK_MODULES([CAIRO], [cairo])])
dnl AS_IF([test "x$enable_phyrex" = "xyes"],[CFLAGS="${CFLAGS} `pkg-config --cflags --libs gtk+-3.0`"])


dnl AM_PATH_GTK_3_0(,,AC_MSG_ERROR(windoe-default needs GTK+-3))

AC_ARG_ENABLE([date],[AS_HELP_STRING([--enable-date],[Compile date])],[date=yes],[date=no])
AM_CONDITIONAL([WANT_DATE], [test "$date" = yes])
if test "$date" = yes; then
   AC_DEFINE([DATE],[1],[DATE tag on])
fi

if test "$phytime" = no; then
if test "$rwrap" = no; then
if test "$test" = no; then
if test "$rf" = no; then
if test "$evolve" = no; then
if test "$phyrex" = no; then
if test "$phyrexsim" = no; then
if test "$date" = no; then
   AC_DEFINE([PHYML],[1],[PHYML tag on])
fi
fi
fi
fi
fi
fi
fi
fi

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT
