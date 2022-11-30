SUFFIX        = f 

DEST	      = .

FFLAGS        = -O3 -ffast-math -floop-block

LIBRARY	      = blas.a

MAKEFILE      = makefile

LINKER        = gfortran

FC            = gfortran

OBJS	      = daxpy.o \
		dcopy.o \
		ddot.o \
		dgemm.o \
		dger.o \
		dscal.o \
		dsymm.o \
		dswap.o \
		dtrsm.o \
		idamax.o \
		lsame.o \
		xerbla.o

SRCS	      = daxpy.f \
		dcopy.f \
		ddot.f \
		dgemm.f \
		dger.f \
		dscal.f \
		dsymm.f \
		dswap.o \
		dtrsm.f \
		idamax.f \
		lsame.f \
		xerbla.f

$(LIBRARY):	$(OBJS); ar -r $(LIBRARY) *.o
		@echo -n "Loading $(LIBRARY) ... "
		@echo "done"
