AC_INIT([randomForestSRC], [3.4.4], [ubk@kogalur.com], [randomForestSRC], [https://github.com/kogalur/randomForestSGT])

AC_PREREQ([2.62])

# Find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi
RBIN="${R_HOME}/bin/R"
LIBS="${PKG_LIBS}"

# Checking for C
# Based on: https://unconj.ca/blog/an-autoconf-primer-for-r-package-authors.html 
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`

AC_LANG(C)
# This the meat of R's m4/openmp.m4
# Based on https://svn.r-project.org/R/trunk/m4/openmp.m4
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
	      for ac_option in -fopenmp -xopenmp -qopenmp \
                         -openmp -mp -omp -qsmp=omp -homp \
                         -fopenmp=libomp \
                         -Popenmp --openmp; do
	        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
	        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
	        AC_LINK_IFELSE([_AC_LANG_OPENMP],
	          [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
	        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
	        if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
	          break
	        fi
	      done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
	      ;; #(
      *)
 	      OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
  fi

# Check for typedefs, structures, and compiler characteristics
AC_CHECK_HEADERS([sys/time.h])

AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(OPENMP_CFLAGS)

# substitute externalized source list
PKG_SOURCES=$(cat src/sources.list)
AC_SUBST(PKG_SOURCES)

AC_CONFIG_FILES([src/Makevars.tmp:src/Makevars.in], [
    if test -f src/Makevars && cmp -s src/Makevars.tmp src/Makevars; then
      AC_MSG_NOTICE([creating src/Makevars])
      AC_MSG_NOTICE([src/Makevars is unchanged])
      rm src/Makevars.tmp
    else
      AC_MSG_NOTICE([creating src/Makevars])
      mv src/Makevars.tmp src/Makevars
    fi
  ]
)

AC_OUTPUT

echo "
 --------------------------------------------------
  Configuration for ${PACKAGE_NAME}
  ================
  
    cflags: ${CFLAGS} 
    libs:   ${PKG_LIBS}

  --------------------------------------------------
"

