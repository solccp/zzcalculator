SRCDIR := src/
OBJDIR := obj/
MODDIR := obj/
BINDIR := bin/

FC := ifort
FFLAGS := -ipo -O3 -funroll-loops -no-prec-div -module ${MODDIR} -pg
#FFLAGS := -g3 -O0 -module ${MODDIR} -trace -check all -debug all

LDFLAGS := -static -lpthread
# -debug all
#LDFLAGS := -static -lmkl_em64t -lguide -lpthread

COMMON_FILES := types_module output
ZHANG_FILES := zhang input functions daughters findZZ polynomial decompose operators print hexagon

COMMON_SRCS := $(addsuffix .f90, ${COMMON_FILES})
COMMON_SRCS := $(addprefix ${SRCDIR}, ${COMMON_SRCS})

COMMON_OBJS := $(addsuffix .o, ${COMMON_FILES})
COMMON_OBJS := $(addprefix ${OBJDIR}, ${COMMON_OBJS})

ZHANG_SRCS := $(addsuffix .f90, ${ZHANG_FILES})
ZHANG_SRCS := $(addprefix ${SRCDIR}, ${ZHANG_SRCS})

ZHANG_OBJS := $(addsuffix .o, ${ZHANG_FILES})
ZHANG_OBJS := $(addprefix ${OBJDIR}, ${ZHANG_OBJS})

SRCS := ${COMMON_SRCS} ${ZHANG_SRCS}

.PHONY: all
all: zhang

.PHONY: zhang
zhang: ${BINDIR}zhang

.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o
	rm -f ${MODDIR}*.mod
	rm .depend

.PHONY: mrproper
mrproper: clean
	rm -f ${BINDIR}*
	rm -f ${SRCDIR}*~
	rm -f *~

${BINDIR}zhang: ${ZHANG_OBJS} ${COMMON_OBJS}
	@echo [LD] $@
	@${FC} ${LDFLAGS} -o $@ $^

${OBJDIR}%.o ${MODDIR}%.mod:${SRCDIR}%.f90
	@echo [FC] $<
	@${FC} ${FFLAGS} -c -o $@ $<

include .depend

.PHONY:	depend
depend .depend: ${SRCS}
	@echo Finding dependencies
	@tools/makedepf90 -b ${OBJDIR} ${SRCS} > .depend

