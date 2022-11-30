SUFFIX        = f 

DEST	      = .

#FFLAGS        = -O3 -ffast-math -floop-block -floop-strip-mine
FFLAGS        = -O3 

LIBRARY	      = nblas.a

MAKEFILE      = makefile

LINKER        = gfortran

FC            = gfortran

OBJS	      = zgesv.o \
		zgetf2.o \
		zgetrf.o \
		zgetrs.o \
		zlaswp.o \
		ztrsm.o \
		xerbla.o \
		dgefa.o \
                dgesl.o \
		ilaenv.o \
		ieeeck.o \
		lsame.o

SRCS	      = zgesv.f \
		zgetf2.f \
		zgetrf.f \
		zgetrs.f \
		zlaswp.f \
		ztrsm.f \
		xerbla.f \
		dgefa.f \
                dgesl.f \
		ilaenv.f \
		ieeeck.f \
		lsame.f

$(LIBRARY):	$(OBJS); ar -r $(LIBRARY) *.o
		@echo -n "Loading $(LIBRARY) ... "
		@echo "done"
