SUFFIX        =.for : f 

DEST	      = .

EXTHDRS	      =

HDRS	      =

!FFLAGS        = -O3 -ffast-math -floop-block
FFLAGS        = -O3 -ffast-math 

!LDFLAGS	      = -O3 -ffast-math -lblas -floop-block
LDFLAGS	      = -O3 -ffast-math -lblas 

LIBS	      = ../../klib/klib.a -lblas

LINKER	      = gfortran-8

FC            = gfortran-8

MAKEFILE      = makefile

OBJS		= modules.o \
		  0kstr.o \
		  alltmp.o \
		  blatts.o \
		  fhndlr.o \
		  gtbess.o \
		  gtneum.o \
		  input.o \
		  jobnam.o \
		  latt3d.o \
		  latt3m.o \
		  layer.o \
		  madl3.o \
		  madl3d.o \
		  mbareh.o \
		  mbares.o \
		  mgaunt.o \
		  primkr.o \
		  primv.o \
		  prnsth.o \
		  prnstr.o \
		  s0lplj.o \
		  s0lplk.o \
		  salpl.o \
		  salplh.o \
		  set3d.o \
		  setcst.o \
		  setflm.o \
		  setplw.o \
		  screen.o \
		  storel.o \
		  storeh.o \
		  stores.o \
		  symsth.o \
		  symstr.o \
		  trmtrx.o \
		  trwats.o \
		  wscell.o

SRCS		= modules.for \
		  0kstr.for \
		  alltmp.for \
		  blatts.for \
		  fhndlr.for \
		  gtbess.for \
		  gtneum.for \
		  input.for \
		  jobnam.for \
		  latt3d.for \
		  latt3m.for \
		  layer.for \
		  madl3.for \
		  madl3d.for \
		  mbareh.for \
		  mbares.for \
		  mgaunt.for \
		  primkr.for \
      	 	  primv.for \
		  prnsth.for \
		  prnstr.for \
		  s0lplj.for \
		  s0lplk.for \
		  salpl.for \
		  salplh.for \
		  set3d.for \
		  setcst.for \
		  setflm.for \
		  setplw.for \
		  screen.for \
		  storel.for \
		  storeh.for \
		  stores.for \
		  symsth.for \
		  symstr.for \
		  trmtrx.for \
		  trwats.for \
		  wscell.for

PRINT	      = pr

PROGRAM	      = kstr

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		@echo -n "Loading $(PROGRAM) ... "
		@$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS)

depend:;	@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@install -s $(PROGRAM) $(DEST)

print:;		@$(PRINT) $(HDRS) $(SRCS)

program:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

.SUFFIXES:      .for $(SUFFIXES)

COMPILE.for=	$(FC) $(FFLAGS) $(TARGET_ARCH) -c
.for.o:
		$(COMPILE.for) $(OUTPUT_OPTION) $<
###
