SUFFIX        = f 

DEST	      = .

LIBRARY       = klib.a

#FFLAGS        = -O3 -ffast-math -floop-block
FFLAGS        = -O3 -fPIC

MAKEFILE      = makefile

LINKER        = gfortran-8

FC            = gfortran-8

OBJS	      = modules.o \
		bessinit.o \
		bessh.o \
		bessl.o \
		besslz.o \
		clbshg.o \
		cross.o \
		dbar.o \
		diffn.o \
		dsifa.o \
		dsidi.o \
		dsisl.o \
		f03ahf.o \
		f03ahz.o \
		f04adf.o \
		f04akf.o \
		f04akz.o \
		fty.o \
		gauntc.o \
		gaunty.o \
		gensim.o \
		grady.o \
		htimer.o \
		htribk.o \
		htridi.o \
		imtql2.o \
		inverc.o \
		kmesh.o \
		p01aaf.o \
		qsort.o \
		realhinit.o \
		realhr.o \
		setclb.o \
		setxcp.o \
		simpn.o \
		volumk.o \
		wagaus.o \
		xch96.o \
		xcpot.o \
		xlda.o

SRCS	      = modules.for \
		bessinit.for \
		bessh.for \
		bessl.for \
		besslz.for \
		clbshg.for \
		cross.for \
		dbar.for \
		diffn.for \
		dsifa.for \
		dsidi.for \
		dsisl.for \
		f03ahf.for \
		f03ahz.for \
		f04adf.for \
		f04akf.for \
		f04akz.for \
		fty.for \
		gauntc.for \
		gaunty.for \
		gensim.for \
		grady.for \
		htimer.for \
		htribk.for \
		htridi.for \
		inverc.for \
		imtql2.for \
		kmesh.for \
		p01aaf.for \
		qsort.for \
		realhinit.for \
		realhr.for \
		setclb.for \
		setxcp.for \
		simpn.for \
		volumk.for \
		wagaus.for \
		xch96.for \
		xcpot.for \
		xlda.for

$(LIBRARY):	$(OBJS); ar -r $(LIBRARY) *.o
		@echo -n "Loading $(LIBRARY) ... "
		@echo "done"
