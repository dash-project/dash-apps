#postfix $(SUBDIRS:%=%-go)
include ./make.defs

MAKEFLAGS = -j 5 -Orecurse

SUBDIRS  = randmat thresh winnow outer product
GOROOT   = baSrcPaper/go
PAPERSUB = expertpar


.PHONY: help printenv all clean all-clean dash dash-clean $(SUBDIRS) 
help:
	@echo "commands:"
	@echo "make            : shows this help message"
	@echo "make help       : shows this help message"
	@echo "make printenv   : prints the environment variables"
	@echo "make dash       : builds all DASH cowichan problems"
	@echo "make dash-clean : cleans all DASH cowichan problems"
	@echo "make problem    : builds the DASH version of the problem. for example: make randmat"
	@echo "make go         : builds all go expertpar cowichan problems"
	@echo "make go-clean   : cleans all go expertpar cowichan problems"
	@echo "make all        : currently builds all of DASH and the expertpar variants"
	@echo "make clean      : currently cleans all of DASH and the expertpar variants"
	@echo "make all-clean  : synonym for clean"
	@echo " "
	@echo "note this makefile is executed with 5 jobs! If you want to change that adjust the '-j' paramter in the MAKEFLAGS."

printenv:
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"
 
 
all: dash go
clean: dash-clean go-clean
all-clean: clean

dash: MAKECMDGOALS = all
dash: $(SUBDIRS)
dash-clean: MAKECMDGOALS = clean
dash-clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

####  <GO>  ####
.PHONY: go go-clean $(SUBDIRS:%=go-%) $(SUBDIRS:%=go-clean-%)


go-clean: $(SUBDIRS:%=go-clean-%)
$(SUBDIRS:%=go-clean-%):
# '-' =>continue on error
	-rm $(GOROOT)/$(@:go-clean-%=%)/expertpar/main

go:$(SUBDIRS:%=go-%)
$(SUBDIRS:%=go-%):
	$(MAKE) -C $(GOROOT)/$(@:go-%=%)/$(PAPERSUB) main
####  </GO>  ####
