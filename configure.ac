AC_INIT([MatrixExtra], 0.1.8)
HAS_LD=`${R_HOME}/bin${R_ARCH_BIN}/Rscript -e "cat(unname(as.integer(capabilities('long.double'))))"`
FLOAT_LIBS=`${R_HOME}/bin${R_ARCH_BIN}/Rscript -e "float:::ldflags(static=Sys.info()[['sysname']] == 'Darwin')"`

AC_LANG(C++)
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[]], [[double *__restrict x;]])],
    [SUPPORTS_RESTRICT="-DSUPPORTS_RESTRICT=1"],
    [SUPPORTS_RESTRICT=""]
)

AC_SUBST(HAS_LD)
AC_SUBST(FLOAT_LIBS)
AC_SUBST(SUPPORTS_RESTRICT)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT 
