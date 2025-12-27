################################################################################
#                                                                              #
#     Make script, Makefile                                                    #
#                                                                              #
#     This script compiles the source programs.                                #
#                                                                              #
#     Author      : Sakakibara Atsushi                                         #
#     Date        : 1998/12/28                                                 #
#     Modification: 1999/01/20, 1999/04/06, 1999/05/20, 1999/06/14,            #
#                   1999/07/05, 1999/11/01, 2000/03/03, 2000/04/18,            #
#                   2000/08/10, 2002/06/17, 2002/07/03, 2002/09/09,            #
#                   2003/05/19, 2004/03/05, 2004/05/31, 2004/09/25,            #
#                   2005/08/05, 2008/10/10, 2009/03/12, 2011/09/22,            #
#                   2013/10/08                                                 #
#                                                                              #
#     Author      : Mayumi Yoshioka                                            #
#     Modification: 2015/03/16                                                 #
#                                                                              #
#     Author      : Masaya Kato                                                #
#     Modification: 2016/04/25                                                 #
#                                                                              #
#     Author      : Satoki Tsujino                                             #
#     Modification: 2017/06/13                                                 #
#                                                                              #
################################################################################

SHELL    = /bin/csh

INCLUDE1 = ${CURDIR}/compile.conf
INCLUDE2 = ${CURDIR}/Src/Make_${TARGET}
INCLUDE3 = ${CURDIR}/Src/Make_object
INCLUDE4 = ${CURDIR}/Src/Make_module
INCLUDE5 = ${CURDIR}/Src/Make_radlib
INCLUDE6 = ${CURDIR}/Src/Make_radcln
INCLUDE7 = ${CURDIR}/Src/Make_spnlib

EXEFILE  = ${TARGET}.exe

LRRTMS   = librrtmg_kgkgb.a
LISPACK  = libisp.a

include ${INCLUDE1}
include ${INCLUDE2}
include ${INCLUDE3}
include ${INCLUDE4}
include ${INCLUDE5}
include ${INCLUDE6}
include ${INCLUDE7}

solver     : ${EXEFILE} ${LRRTMS}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

radlib    : ${LRRTMS}
	@${TOUCH} $@
	@echo -n 'Last library object was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

spnlib    : ${LISPACK}
	@${TOUCH} $@
	@echo -n 'Last library object was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

gridata    : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@
asldata    : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

radata     : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

terrain    : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

surface    : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

unite      : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

rstruct    : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

check      : ${EXEFILE}
	@${TOUCH} $@
	@echo -n 'Last executable was generated at ' >> $@
	@echo `${DATE} "+%Y/%m/%d %T"`'.' >> $@

clean      :
	${RM} -f ${EXES} ${MODS} ${ADDMODS} ${OBJS} ${ADDOBJS}
	${RM} -f ${TMPSRC} ${LISOBJS} *.s
	${RM} -f ${RADOBJS} ${RADMODS} ${LRRTMS} ${LISPACK}

crscln     :
	${RM} -f ${EXES} ${MODS} ${ADDMODS} ${OBJS} ${ADDOBJS}

radcln     :
	${RM} -f ${RADOBJS} ${RADMODS} ${LRRTMS}

spncln     :
	${RM} -f ${LISOBJS} ${TMPSRC} ${LISPACK}

ifeq ($(TARGET), solver)
${EXEFILE} : ${OBJS} ${LISOBJS} ${ADDOBJS} ${MAINOBJ} ${LRRTMS} ${LISPACK}
	${FC} -o $@ ${MAINOBJ} ${OBJS} ${ADDOBJS} ${LDFLAGS} ${LOCALLIBS}
else
${EXEFILE} : ${OBJS} ${ADDOBJS} ${MAINOBJ}
	${FC} -o $@ ${MAINOBJ} ${OBJS} ${ADDOBJS} ${LDFLAGS}
endif

${LRRTMS} : ${LRROBJS}
	${AR} -r ${LRRTMS} ${LRROBJS}

${LISPACK} : ${TMPSRC} ${LISOBJS}
	$(AR) $(ARFLAGS) $(STATICLIB) ${LISOBJS}
	$(RANLIB) $(STATICLIB)

.SUFFIXES  :

.SUFFIXES  : .o .f90 .c .f .s .F

.f90.o     :
	${FC} -c ${FFLAGS} $<

.c.o       :
	${CC} -c ${CFLAGS} $<

.F.o :
	$(CPP) -C -traditional $(CPPFLAGS) $< | sed '/^#/d' > $*.f
	$(FC) $(FFLAGS) -c $*.f
	
.f.o :
	$(FC) $(FFLAGS) -c $<
	
.s.o :
	$(AS) $< -o $@
	
.F.f :
	$(CPP) -C -traditional $(CPPFLAGS) $< | sed '/^#/d' > $*.f
