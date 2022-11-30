SUFFIX        = .f

DEST	      = .

EXTHDRS	      =

HDRS	      =

FFLAGS        = -O3 -ffast-math -floop-block

LDFLAGS	      = -O3 -ffast-math -lblas -floop-block

LIBS	      = ../klib/klib.a -lblas

LINKER	      = gfortran

FC            = gfortran

MAKEFILE      = makefile

OBJS		= modules.o \
		  0bmdl.o \
		  fhndlr.o  \
		  gaunt.o \
		  gamfc.o \
		  jobnam.o \
		  latt3m.o \
		  madl3m.o \
		  primv.o \
		  primkr.o \
		  set3d.o \
		  setcst.o \
		  setup.o \
		  smtrx.o \
		  ylmrk.o

SRCS		= modules.for \
		  0bmdl.for \
		  fhndlr.for  \
		  gamfc.for \
		  gaunt.for \
		  jobnam.for \
		  latt3m.for \
		  madl3m.for \
		  primv.for \
		  primkr.for \
		  set3d.for \
		  setcst.for \
		  setup.for \
		  smtrx.for \
		  ylmrk.for

PRINT	      = pr

PROGRAM	      = bmdl

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
