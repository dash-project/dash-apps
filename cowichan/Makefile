#suffix $(SUBDIRS:%=%-go)
#@echo $(@:%-go=%) #get rid of suffix
include ./make.defs

SHELL:=/bin/bash

GO_ROOT   = baSrcPaper/go
CHPL_ROOT = baSrcPaper/chapel
CILK_ROOT = baSrcPaper/cilk
TBB_ROOT  = baSrcPaper/tbb
SUBDIRS   = randmat thresh winnow outer product


####################         <Info>          ####################
.PHONY: help printenv 
#### ------------------------------------------------------------
help:
	@echo "commands:"
	@echo "make            : shows this help message"
	@echo "make help       : shows this help message"
	@echo "make printenv   : prints the environment variables"
	@echo "make dash       : builds all DASH   cowichan problems"
	@echo "make dash-clean : cleans all DASH   cowichan problems"
	@echo "make go         : builds all go     expertpar cowichan problems"
	@echo "make go-clean   : cleans all go     expertpar cowichan problems"
	@echo "make chpl       : builds all chapel expertpar cowichan problems"
	@echo "make chpl-clean : cleans all chapel expertpar cowichan problems"
	@echo "make cilk       : builds all cilk   expertpar cowichan problems"
	@echo "make cilk-clean : cleans all cilk   expertpar cowichan problems"
	@echo "make tbb        : builds all tbb    expertpar cowichan problems"
	@echo "make tbb-clean  : cleans all tbb    expertpar cowichan problems"
	@echo "***************************************************************"
	@echo "make all        : builds all of DASH and the expertpar variants"
	@echo "make clean      : cleans all of DASH and the expertpar variants"
	@echo "make all-clean  : synonym for clean"

printenv:
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"
####################        </Info>          ####################  
####################  <Top Level Targets>     ###################
.PHONY: all clean all-clean
#### ------------------------------------------------------------
all: dash-gnu go chpl cilk tbb

clean: dash-clean go-clean chpl-clean cilk-clean tbb-clean
all-clean: clean
####################  </Top Level Targets>  #####################
####################       <DASH>       #########################
.PHONY: dash dash-clean $(SUBDIRS)
#### ------------------------------------------------------------
dash-gnu:
ifeq (icc, $(findstring icc, $(CC)))
	@module swap intel gnu; \
	$(MAKE) dash
else
	$(MAKE) dash
endif

dash-intel:
ifeq (gcc, $(CC))
	@module swap gnu intel; \
	$(MAKE) dash
else
	$(MAKE) dash
endif

dash: MAKECMDGOALS = all
dash: $(SUBDIRS)

dash-clean: MAKECMDGOALS = clean
dash-clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)
####################    </DASH>    ##############################
####################     <GO>      ##############################
.PHONY: go go-clean $(SUBDIRS:%=go-%) $(SUBDIRS:%=go-clean-%)
#### ------------------------------------------------------------
go-clean: $(SUBDIRS:%=go-clean-%)
$(SUBDIRS:%=go-clean-%):
	-rm $(GO_ROOT)/$(@:go-clean-%=%)/expertpar/main
	#-rm $(GO_ROOT)/$(@:go-clean-%=%)/expertseq/main

go:$(SUBDIRS:%=go-%)
$(SUBDIRS:%=go-%):
	$(MAKE) -C $(GO_ROOT)/$(@:go-%=%)/expertpar main
	#$(MAKE) -C $(GO_ROOT)/$(@:go-%=%)/expertseq main
####################     </GO>     ##############################
####################    <Chapel>   ##############################
.PHONY: chpl chpl-clean $(SUBDIRS:%=chpl-%) $(SUBDIRS:%=chpl-clean-%) chpl-test
#### ------------------------------------------------------------
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
####################    </Chapel>  ##############################
####################     <Cilk>    ##############################
.PHONY: cilk cilk_ cilk-clean $(SUBDIRS:%=cilk-%) $(SUBDIRS:%=cilk-clean-%)
#### ------------------------------------------------------------
cilk:
ifeq (gcc, $(CC))
	@module swap gnu intel; \
	$(MAKE) cilk_
else
	$(MAKE) cilk_
endif

cilk_: $(SUBDIRS:%=cilk-%)

$(SUBDIRS:%=cilk-%):
	$(MAKE) -C $(CILK_ROOT)/$(@:cilk-%=%)/expertpar main

cilk-clean: $(SUBDIRS:%=cilk-clean-%)
$(SUBDIRS:%=cilk-clean-%):
	-rm $(CILK_ROOT)/$(@:cilk-clean-%=%)/expertpar/main
####################    </Cilk>    ##############################
####################     <TBB>    ##############################
.PHONY: tbb tbb_ tbb-clean $(SUBDIRS:%=tbb-%) $(SUBDIRS:%=tbb-clean-%)
#### ------------------------------------------------------------
tbb:
ifeq (gcc, $(CC))
	@module swap gnu intel; \
	$(MAKE) tbb_
else
	$(MAKE) tbb_
endif

tbb_: $(SUBDIRS:%=tbb-%)

$(SUBDIRS:%=tbb-%):
	$(MAKE) -C $(TBB_ROOT)/$(@:tbb-%=%)/expertpar main

tbb-clean: $(SUBDIRS:%=tbb-clean-%)
$(SUBDIRS:%=tbb-clean-%):
	-rm $(TBB_ROOT)/$(@:tbb-clean-%=%)/expertpar/main
####################    </TBB>    ##############################
