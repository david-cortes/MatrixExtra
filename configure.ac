AC_INIT([MatrixExtra], 0.1.6)
HAS_LD=`${R_HOME}/bin${R_ARCH_BIN}/Rscript -e "cat(unname(as.integer(capabilities('long.double'))))"`
FLOAT_LIBS=`${R_HOME}/bin${R_ARCH_BIN}/Rscript -e "float:::ldflags(static=Sys.info()[['sysname']] == 'Darwin')"`

AC_SUBST(HAS_LD)
AC_SUBST(FLOAT_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT 
