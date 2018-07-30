#-----------------------------------------------------------------------
# File    : makefile
# Contents: build sequoia program (on Unix systems)
# Author  : Christian Borgelt
# History : 2010.08.06 file created from relim makefile
#           2010.08.22 module escape added (for module tabread)
#           2010.10.08 changed standard from -ansi to -std=c99
#           2013.03.20 extended the requested warnings in CFBASE
#           2013.10.19 modules tabread and patspec added
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../sequoia/src
UTILDIR  = ../../util/src
TRACTDIR = ../../tract/src

CC       = gcc -std=c99
# CC       = g++
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic -c $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(TRACTDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
           $(UTILDIR)/symtab.h   $(UTILDIR)/error.h    \
           $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
           $(TRACTDIR)/tract.h   $(TRACTDIR)/patspec.h \
           $(TRACTDIR)/report.h
OBJS     = $(UTILDIR)/arrays.o   $(UTILDIR)/idmap.o    \
           $(UTILDIR)/escape.o   $(UTILDIR)/tabread.o  \
           $(UTILDIR)/tabwrite.o $(UTILDIR)/scform.o   \
           $(TRACTDIR)/taread.o  $(TRACTDIR)/patspec.o \
           $(TRACTDIR)/report.o  sequoia.o $(ADDOBJS)
PRGS     = sequoia

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

sequoia:      $(OBJS) makefile
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
sequoia.o:    $(HDRS)
sequoia.o:    sequoia.c makefile
	$(CC) $(CFLAGS) $(INCS) -DSOIA_MAIN sequoia.c -o $@

sequoia.d:    sequoia.c
	$(CC) -MM $(CFLAGS) $(INCS) -DSOIA_MAIN sequoia.c > sequoia.d

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/idmap.o:
	cd $(UTILDIR);  $(MAKE) idmap.o    ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/escape.o:
	cd $(UTILDIR);  $(MAKE) escape.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabread.o:
	cd $(UTILDIR);  $(MAKE) tabread.o  ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabwrite.o:
	cd $(UTILDIR);  $(MAKE) tabwrite.o ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/scform.o:
	cd $(UTILDIR);  $(MAKE) scform.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/taread.o:
	cd $(TRACTDIR); $(MAKE) taread.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/patspec.o:
	cd $(TRACTDIR); $(MAKE) patspec.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/report.o:
	cd $(TRACTDIR); $(MAKE) report.o   ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f sequoia.zip sequoia.tar.gz; \
        zip -rq sequoia.zip sequoia/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz sequoia.tar.gz sequoia/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) /usr/local/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.d *.o *~ *.flc core $(PRGS)

clean:
	$(MAKE) localclean
	cd $(UTILDIR);  $(MAKE) clean
	cd $(TRACTDIR); $(MAKE) localclean
