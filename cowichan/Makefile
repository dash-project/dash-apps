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
	@echo "make go         : builds all go     expertpar cowichan problems"
	@echo "make chpl       : builds all chapel expertpar cowichan problems"
	@echo "make cilk       : builds all cilk   expertpar cowichan problems"
	@echo "make tbb        : builds all tbb    expertpar cowichan problems"
	@echo "make dash-clean : cleans all DASH   cowichan problems"
	@echo "make   go-clean : cleans all go     expertpar cowichan problems"
	@echo "make chpl-clean : cleans all chapel expertpar cowichan problems"
	@echo "make cilk-clean : cleans all cilk   expertpar cowichan problems"
	@echo "make  tbb-clean : cleans all tbb    expertpar cowichan problems"
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
####################   <Top Level Targets>   ####################
.PHONY: all clean all-clean
#### ------------------------------------------------------------
all: dash-gnu go chpl cilk tbb

clean: dash-clean go-clean chpl-clean cilk-clean tbb-clean test-clean
all-clean: clean
####################  </Top Level Targets>  #####################
####################       <DASH>           #####################
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
chpl:
ifeq (icc, $(findstring icc, $(CC)))
	@module swap intel gnu; \
	$(MAKE) chpl_
else
	$(MAKE) chpl_
endif

chpl_: $(SUBDIRS:%=chpl-%)

$(SUBDIRS:%=chpl-%):
	$(MAKE) -C $(CHPL_ROOT)/$(@:chpl-%=%)/expertpar main

chpl-clean: $(SUBDIRS:%=chpl-clean-%)
$(SUBDIRS:%=chpl-clean-%):
	-rm $(CHPL_ROOT)/$(@:chpl-clean-%=%)/expertpar/main
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
####################     <TBB>     ##############################
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
####################    </TBB>     ##############################
################   <Validation Test>   ##########################
.PHONY: test-all

test_results:; mkdir -p test_results

test: test-all
test-all: chpl-test go-test cilk-test tbb-test dash-test | test_results
	@echo 
	@echo Comparing outputs...
	@$(MAKE) --no-print-directory compare
# $(CHPL_ROOT)/thresh/expertpar/main  < test_results/chpl_randmat

.PHONY: dash-test dash-test-randmat
dash-test: | test_results
ifeq (icc, $(findstring icc, $(CC)))
	@module swap intel gnu &> /dev/null; \
	$(MAKE) --no-print-directory dash-test-randmat
	@echo DASH.....finished
else
	@$(MAKE) --no-print-directory dash-test-randmat
	@echo DASH.....finished
endif
dash-test-randmat: | test_results; @echo "10 9 8" | ./randmat/randmat > test_results/dash_randmat_out

.PHONY: chpl-test chpl-test-randmat
chpl-test:| test_results
ifeq (icc, $(findstring icc, $(CC)))
	@module swap intel gnu &> /dev/null; \
	$(MAKE) --no-print-directory chpl-test-randmat
	@echo Chapel...finished
else
	@$(MAKE) --no-print-directory chpl-test-randmat
	@echo Chapel...finished
endif
chpl-test-randmat: | test_results; @echo "10 9 8" | $(CHPL_ROOT)/randmat/expertpar/main > test_results/chpl_randmat_out

.PHONY: go-test go-test-randmat
go-test: go-test-randmat | test_results
	@echo Go.......finished
go-test-randmat:|test_results; @echo "10 9 8" | $(GO_ROOT)/randmat/expertpar/main > test_results/go_randmat_out

.PHONY: cilk-test cilk-test-randmat
cilk-test: |test_results
ifeq (gcc,$(CC))
	@module swap gnu intel &> /dev/null; \
	$(MAKE) --no-print-directory cilk-test-randmat
	@echo Cilk.....finished
else
	@$(MAKE) --no-print-directory cilk-test-randmat
	@echo Cilk.....finished
endif
cilk-test-randmat:|test_results; @echo "10 9 8" | $(CILK_ROOT)/randmat/expertpar/main > test_results/cilk_randmat_out

.PHONY: tbb-test tbb-test-randmat
tbb-test:| test_results
ifeq (gcc,$(CC))
	@module swap gnu intel &> /dev/null; \
	$(MAKE) --no-print-directory tbb-test-randmat
	@echo TBB......finished
else
	@$(MAKE) --no-print-directory tbb-test-randmat
	@echo TBB......finished
endif
tbb-test-randmat:|test_results; @echo "10 9 8" | $(TBB_ROOT)/randmat/expertpar/main > test_results/tbb_randmat_out

.PHONY: compare cmp-randmat
compare: cmp-randmat
cmp-randmat: RANDMAT_RESULT = $(if\
		$(and \
			$(findstring identical, $(shell diff -sw test_results/dash_randmat_out test_results/chpl_randmat_out)),\
			$(findstring identical, $(shell diff -sw test_results/dash_randmat_out test_results/go_randmat_out)),\
			$(findstring identical, $(shell diff -sw test_results/dash_randmat_out test_results/cilk_randmat_out)),\
			$(findstring identical, $(shell diff -sw test_results/dash_randmat_out test_results/tbb_randmat_out))\
		),\
	SUCCESS, FAILURE)
cmp-randmat:;	@echo Randmat Validation: $(RANDMAT_RESULT)
	
.PHONY: test-clean
test-clean:
	-rm -r test_results
################  </Validation Test>   ##########################
