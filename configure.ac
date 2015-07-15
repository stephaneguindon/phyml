#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([PhyML],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[s.guindon@auckland.ac.nz])
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

AC_CONFIG_SRCDIR([src/simu.c],[doc/phyml-manual.tex])
AC_CONFIG_HEADERS([config.h])

AC_DEFINE([UNIX],[1],[Unix tag on])
AC_DEFINE([DEBUG],[1],[Debug tag on])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

AC_PROG_CC([gcc])


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

# Thanks for Charles Plessy from Debian
case "${host}" in
*i386*|*amd64*|*x86_64*)
MSSE=-msse
AC_MSG_NOTICE("On arch ${host} use MSSE = ${MSSE}")
;;
*)
MSSE=
AC_MSG_NOTICE("Do not use MSSE = ${MSSE} because not available on arch ${host}")
;;
esac

dnl CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -Wall ${ARCH_flag}"
dnl CFLAGS="-O3 -Wfloat-equal -fomit-frame-pointer -funroll-loops"
dnl CFLAGS="-O3 -Wfloat-equal"
dnl CFLAGS="-g"

LDFLAGS="${ARGC_flag}" 

LT_INIT
AC_PROG_LIBTOOL
dnl Add option for compiling with debugging informations

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Remove optimization options and add debug informations.])])
AS_IF([test "x$enable_debug" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g -Wuninitialized"],
      [CFLAGS="-Wall -O2 ${MSSE} -fomit-frame-pointer -funroll-loops -Wempty-body -Wuninitialized ${ARCH_flag}"])

AC_ARG_ENABLE([gprof],
              [AS_HELP_STRING([--enable-gprof],
                              [Remove optimization options and add profiling informations.])])
AS_IF([test "x$enable_gprof" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g -pg"],
      [AS_IF([test "x$enable_debug" = "xno"],[CFLAGS="-Wall -O2 -msse -fomit-frame-pointer -funroll-loops ${ARCH_flag}"])])

AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
                              [Compile with mpicc instead of gcc.])])
AS_IF([test "x$enable_mpi" = "xyes"],[CC="mpicc"])
AS_IF([test "x$enable_mpi" = "xyes"],AC_DEFINE([MPI],[1],[MPI tag on]))
AM_CONDITIONAL([WANT_MPI], [test "x$enable_mpi" = "xyes"])


AC_ARG_ENABLE([win],
              [AS_HELP_STRING([--enable-win],
                              [Compile with mingw instead of gcc.])])
AS_IF([test "x$enable_win" = "xyes"],[CC="i686-w64-mingw32-gcc"])


AC_ARG_ENABLE([beagle], [AS_HELP_STRING([--enable-beagle], [Compute likelihoods using BEAGLE library.])], [beagle=yes],[beagle=no])
AS_IF([test "x$enable_beagle" = "xyes"],[PKG_CHECK_MODULES(BEAGLE, hmsbeagle-1)])
AS_IF([test "x$enable_beagle" = "xyes"],
      [CFLAGS="${CFLAGS} ${BEAGLE_CFLAGS} ${BEAGLE_LIBS} -lstdc++ ${ARCH_flag}"],
      [CFLAGS=${CFLAGS}])
AS_IF([test "x$enable_beagle" = "xyes"], AC_DEFINE([BEAGLE],[1],[BEAGLE tag on]))
AM_CONDITIONAL([WANT_BEAGLE], [test "$enable_beagle" = yes])


AC_ARG_ENABLE([phytime],[AS_HELP_STRING([--enable-phytime],[Compile PhyTime])],[phytime=yes],[phytime=no])
AM_CONDITIONAL([WANT_PHYTIME], [test "$phytime" = yes])
if test "$phytime" = yes; then
   AC_DEFINE([PHYTIME],[1],[PHYTIME tag on])
fi

AC_ARG_ENABLE([tiporder],[AS_HELP_STRING([--enable-tiporder],[Compile tiporder])],[tiporder=yes],[tiporder=no])
AM_CONDITIONAL([WANT_TIPORDER], [test "$tiporder" = yes])
if test "$tiporder" = yes; then
   AC_DEFINE([TIPORDER],[1],[TIPORDER tag on])
fi

AC_ARG_ENABLE([part],[AS_HELP_STRING([--enable-part],[Compile Part])],[part=yes],[part=no])
AM_CONDITIONAL([WANT_PART], [test "$part" = yes])
if test "$part" = yes; then
   AC_DEFINE([PART],[1],[PART tag on])
fi

AC_ARG_ENABLE([rwrap],[AS_HELP_STRING([--enable-rwrap],[Compile Rwrap])],[rwrap=yes],[rwrap=no])
AM_CONDITIONAL([WANT_RWRAP], [test "$rwrap" = yes])
if test "$rwrap" = yes; then
   AC_DEFINE([RWRAP],[1],[RWRAP tag on])
fi

AC_ARG_ENABLE([phycont],[AS_HELP_STRING([--enable-phycont],[Compile PhyCont])],[phycont=yes],[phycont=no])
AM_CONDITIONAL([WANT_PHYCONT], [test "$phycont" = yes])
if test "$phycont" = yes; then
   AC_DEFINE([PHYCONT],[1],[PHYCONT tag on])
fi

AC_ARG_ENABLE([m4],[AS_HELP_STRING([--enable-m4],[Compile M4])],[m4=yes],[m4=no])
AM_CONDITIONAL([WANT_M4], [test "$m4" = yes])
if test "$m4" = yes; then
   AC_DEFINE([M4],[1],[m4 tag on])
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

AC_ARG_ENABLE([invitee],[AS_HELP_STRING([--enable-invitee],[Compile invitee])],[invitee=yes],[invitee=no])
AM_CONDITIONAL([WANT_INVITEE], [test "$invitee" = yes])
if test "$invitee" = yes; then
   AC_DEFINE([INVITEE],[1],[INVITEE tag on])
fi

AC_ARG_ENABLE([geo],[AS_HELP_STRING([--enable-geo],[Compile geo])],[geo=yes],[geo=no])
AM_CONDITIONAL([WANT_GEO], [test "$geo" = yes])
if test "$geo" = yes; then
   AC_DEFINE([GEO],[1],[GEO tag on])
fi

AC_ARG_ENABLE([checkpoint],[AS_HELP_STRING([--enable-checkpoint],[Compile checkpoint])],[checkpoint=yes],[checkpoint=no])
AM_CONDITIONAL([WANT_CHECKPOINT], [test "$checkpoint" = yes])
if test "$checkpoint" = yes; then
   AC_DEFINE([CHECKPOINT],[1],[CHECKPOINT tag on])
fi


AC_ARG_ENABLE([migrep],[AS_HELP_STRING([--enable-migrep],[Compile migrep])],[migrep=yes],[migrep=no])
AM_CONDITIONAL([WANT_MIGREP], [test "$migrep" = yes])
if test "$migrep" = yes; then
   AC_DEFINE([MIGREP],[1],[MIGREP tag on])
fi
dnl AS_IF([test "x$enable_migrep" = "xyes"],[PKG_CHECK_MODULES([GTK], [gtk+-3.0])])
dnl AS_IF([test "x$enable_migrep" = "xyes"],[PKG_CHECK_MODULES([CAIRO], [cairo])])
dnl AS_IF([test "x$enable_migrep" = "xyes"],[CFLAGS="${CFLAGS} `pkg-config --cflags --libs gtk+-3.0` ${ARCH_flag}"])


dnl AM_PATH_GTK_3_0(,,AC_MSG_ERROR(windoe-default needs GTK+-3))

if test "$phytime" = no; then
if test "$tiporder" = no; then
if test "$part" = no; then
if test "$rwrap" = no; then
if test "$phycont" = no; then
if test "$m4" = no; then
if test "$test" = no; then
if test "$invitee" = no; then
if test "$geo" = no; then
if test "$evolve" = no; then
if test "$checkpoint" = no; then
if test "$migrep" = no; then
   AC_DEFINE([PHYML],[1],[PHYML tag on])
fi
fi
fi
fi
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
