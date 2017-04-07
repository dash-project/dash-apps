#suffix $(SUBDIRS:%=%-go)
#@echo $(@:%-go=%) #get rid of suffix
include ./make.defs

GO_ROOT   = baSrcPaper/go
CHPL_ROOT = baSrcPaper/chapel
CILK_ROOT = baSrcPaper/cilk
SUBDIRS   = randmat thresh winnow outer product

SHELL:=/bin/bash

.PHONY: help printenv 
help:
	@echo "commands:"
	@echo "make            : shows this help message"
	@echo "make help       : shows this help message"
	@echo "make printenv   : prints the environment variables"
	@echo "make dash       : builds all DASH cowichan problems"
	@echo "make dash-clean : cleans all DASH cowichan problems"
	@echo "make problem    : builds the DASH version of the problem. for example: make randmat"
	@echo "make go         : builds all go     expertpar/seq cowichan problems"
	@echo "make go-clean   : cleans all go     expertpar/seq cowichan problems"
	@echo "make chpl       : builds all chapel expertpar/seq cowichan problems"
	@echo "make chpl-clean : cleans all chapel expertpar/seq cowichan problems"
	@echo "make all        : currently builds all of DASH and the expertpar/seq variants"
	@echo "make clean      : currently cleans all of DASH and the expertpar/seq variants"
	@echo "make all-clean  : synonym for clean"

printenv:
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.PHONY: all clean all-clean dash dash-clean init-modules swap-gnu swap-intel $(SUBDIRS)

init-modules:
	source /etc/profile.d/modules.sh

swap-gnu: init-modules
	module swap intel gnu

swap-intel: init-modules
	module swap gnu intel

all: swap-gnu dash go chpl #intel
all-clean: clean
clean: dash-clean go-clean chpl-clean cilk-clean


.PHONY: intel gnu testB testA
testA:
	@echo $(PATH)

testB:
	@echo "testB hier"
	@echo $(CC)

gnu:
ifneq (gcc, $(CC)) 
	./swap.sh
	@echo "gnu hier!"
	@echo $(CC)
#module swap intel gnu
else
	@echo bla
	@echo $(CC)
endif

intel: gnu
#	$(MAKE)


####  <DASH>  ####
dash: MAKECMDGOALS = all
dash: $(SUBDIRS)

dash-gnu: swap-gnu dash

dash-intel: swap-intel dash

dash-clean: MAKECMDGOALS = clean
dash-clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)
####  </DASH>  ####

####  <GO>  ####
.PHONY: go go-clean $(SUBDIRS:%=go-%) $(SUBDIRS:%=go-clean-%)

go-clean: $(SUBDIRS:%=go-clean-%)
$(SUBDIRS:%=go-clean-%):
	-rm $(GO_ROOT)/$(@:go-clean-%=%)/expertpar/main
	-rm $(GO_ROOT)/$(@:go-clean-%=%)/expertseq/main

go:$(SUBDIRS:%=go-%)
$(SUBDIRS:%=go-%):
	$(MAKE) -C $(GO_ROOT)/$(@:go-%=%)/expertpar main
	$(MAKE) -C $(GO_ROOT)/$(@:go-%=%)/expertseq main
####  </GO>  ####


#### <Chapel> ####
.PHONY: chpl chpl-clean $(SUBDIRS:%=chpl-%) $(SUBDIRS:%=chpl-clean-%) chpl-test
chpl: $(SUBDIRS:%=chpl-%)

$(SUBDIRS:%=chpl-%):
	$(MAKE) -C $(CHPL_ROOT)/$(@:chpl-%=%)/expertpar main

chpl-clean: $(SUBDIRS:%=chpl-clean-%)
$(SUBDIRS:%=chpl-clean-%):
	-rm $(CHPL_ROOT)/$(@:chpl-clean-%=%)/expertpar/main

chpl-test: #chpl-randmat chpl-thresh
	mkdir -p test_results
	@echo "10 9 8" | $(CHPL_ROOT)/randmat/expertpar/main > test_results/chpl_randmat
# $(CHPL_ROOT)/thresh/expertpar/main  < test_results/chpl_randmat
#### </Chapel> ####


#### <Cilk> ####
.PHONY: cilk cilk-clean $(SUBDIRS:%=cilk-%) $(SUBDIRS:%=cilk-clean-%)
cilk: $(SUBDIRS:%=cilk-%)

$(SUBDIRS:%=cilk-%):
	$(MAKE) -C $(CILK_ROOT)/$(@:cilk-%=%)/expertpar main

cilk-clean: $(SUBDIRS:%=cilk-clean-%)
$(SUBDIRS:%=cilk-clean-%):
	-rm $(CILK_ROOT)/$(@:cilk-clean-%=%)/expertpar/main
#### </Cilk> ####
