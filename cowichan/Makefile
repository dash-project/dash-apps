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

clean: dash-clean go-clean chpl-clean cilk-clean tbb-clean run-clean
all-clean: clean
####################  </Top Level Targets>  #####################
####################       <DASH>           #####################
.PHONY: dash-gnu dash-intel dash dash-clean $(SUBDIRS)
#### ------------------------------------------------------------
dash-gnu:
	@if [[ $$CC == *icc ]] ; then module swap intel gnu; fi;\
	$(MAKE) dash

dash-intel:
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) dash

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

go:$(SUBDIRS:%=go-%)
$(SUBDIRS:%=go-%):
	$(MAKE) -C $(GO_ROOT)/$(@:go-%=%)/expertpar main
####################     </GO>     ##############################
####################    <Chapel>   ##############################
.PHONY: chpl chpl_ chpl-clean $(SUBDIRS:%=chpl-%) $(SUBDIRS:%=chpl-clean-%)
#### ------------------------------------------------------------
chpl:
	@if [[ $$CC == *icc ]] ; then module swap intel gnu; fi;\
	$(MAKE) chpl_

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
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) cilk_

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
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) tbb_

tbb_: $(SUBDIRS:%=tbb-%)

$(SUBDIRS:%=tbb-%):
	$(MAKE) -C $(TBB_ROOT)/$(@:tbb-%=%)/expertpar main

tbb-clean: $(SUBDIRS:%=tbb-clean-%)
$(SUBDIRS:%=tbb-clean-%):
	-rm $(TBB_ROOT)/$(@:tbb-clean-%=%)/expertpar/main
####################    </TBB>     ##############################
################   <Validation Test>   ##########################
run_results:; mkdir -p run_results

.PHONY: run run-all
run: run-all
run-all: chpl-run go-run cilk-run tbb-run dash-run | run_results

.PHONY: dash-run dash-run_ dash-run-randmat
dash-run: | run_results
	@if [[ $$CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;\
	echo DASH ..... run with: rows:$(rows), cols:$(cols), seed:$(seed);\
	$(MAKE) --no-print-directory dash-run_ ;\
	echo DASH ..... finished

dash-run_:
	@echo $(rows) $(cols) > run_results/dash_randmat_out
	@echo $(rows) $(cols) $(seed) | ./randmat/randmat >> run_results/dash_randmat_out

	
.PHONY: chpl-run chpl-run_ chpl-run-randmat
chpl-run:| run_results
	@if [[ $$CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;\
	echo Chapel ... run with: rows:$(rows), cols:$(cols), seed:$(seed);\
	$(MAKE) --no-print-directory chpl-run_ ;\
	echo Chapel ... finished

chpl-run_:
	@echo $(rows) $(cols) > run_results/chpl_randmat_out
	@echo $(rows) $(cols) $(seed) | $(CHPL_ROOT)/randmat/expertpar/main >> run_results/chpl_randmat_out

	
.PHONY: go-run go-run_ go-run-randmat
go-run:| run_results
	@echo Go ....... run with: rows:$(rows), cols:$(cols), seed:$(seed)
	@$(MAKE) --no-print-directory go-run_
	@echo Go ....... finished

go-run_: 
	@echo $(rows) $(cols) > run_results/go_randmat_out
	@echo $(rows) $(cols) $(seed) | $(GO_ROOT)/randmat/expertpar/main >> run_results/go_randmat_out

	
.PHONY: cilk-run cilk-run_ cilk-run-randmat
cilk-run: |run_results
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;\
	echo Cilk ..... run with: rows:$(rows), cols:$(cols), seed:$(seed);\
	$(MAKE) --no-print-directory cilk-run_ ;\
	echo Cilk ..... finished

cilk-run_:
	@echo $(rows) $(cols) > run_results/cilk_randmat_out
	@echo $(rows) $(cols) $(seed) | $(CILK_ROOT)/randmat/expertpar/main >> run_results/cilk_randmat_out


.PHONY: tbb-run tbb-run_ tbb-run-randmat
tbb-run:| run_results
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;\
	echo TBB ...... run with: rows:$(rows), cols:$(cols), seed:$(seed);\
	$(MAKE) --no-print-directory tbb-run_ ;\
	echo TBB ...... finished

tbb-run_:
	@echo $(rows) $(cols) > run_results/tbb_randmat_out
	@echo $(rows) $(cols) $(seed) | $(TBB_ROOT)/randmat/expertpar/main >> run_results/tbb_randmat_out
# @echo $(rows) $(cols) $(thresh) | $(TBB_ROOT)/thresh/expertpar/main > run_results/tbb_thresh_out


RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

SUCCESS:=$(GREEN)SUCCESS$(NC)
FAILURE:=$(RED)FAILURE$(NC)

.PHONY: test
test:
	@echo How many iterations should be made? \(Enter a count and press \[Enter\]\):; read i ;\
	echo ---------- ;\
	for (( succ=0 , fail=0 , subsucc=0 ; i-- ; subsucc=0 )) ; do \
	export rows=$$(( 1 + RANDOM % 200 )) cols=$$(( 1 + RANDOM % 200 )) seed=$$(( 1 + RANDOM % 60 )) thresh=$$(( 10 + RANDOM % 90 ));\
		$(MAKE) --no-print-directory run-all;\
		echo comparing ... ;\
		if [[ $$(diff -sw run_results/dash_randmat_out run_results/chpl_randmat_out) == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)Chapel randmat output differs from DASH!$(NC); fi;\
		if [[ $$(diff -sw run_results/dash_randmat_out run_results/go_randmat_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)GO randmat output differs from DASH!$(NC); fi;\
		if [[ $$(diff -sw run_results/dash_randmat_out run_results/cilk_randmat_out) == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)Cilk randmat output differs from DASH!$(NC); fi;\
		if [[ $$(diff -sw run_results/dash_randmat_out run_results/tbb_randmat_out)  == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)TBB randmat output differs from DASH!$(NC); fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)randmat output comparison without differences.$(NC);;\
			*) (( ++fail )) ;;\
		esac;\
	done ;\
	echo -e $(GREEN)successes$(NC):$$succ $(RED)failures$(NC):$$fail ;\
	
.PHONY: run-clean
run-clean:
	-rm -r run_results


.PHONY: testA testB
testA:
	@if [[ $$(diff -sw run_results/dash_randmat_out run_results/chpl_randmat_out) == *are\ identical* ]] ;then echo true; else echo false; fi;
	@if [[ $$(diff -sw run_results/dash_randmat_out run_results/tbb_randmat_out) == *are\ identical* ]] ;then echo true; else echo false; fi;


testB:
	@echo $$((2#1001)) $$((b=~2#1001)) $$b
	
################  </Validation Test>   ##########################
