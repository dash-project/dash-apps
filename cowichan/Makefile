#suffix $(SUBDIRS:%=%-go)
#@echo $(@:%-go=%) #get rid of suffix
include ./make.defs

SHELL:=/bin/bash

GO_ROOT   = baSrcPaper/go
CHPL_ROOT = baSrcPaper/chapel
CILK_ROOT = baSrcPaper/cilk
TBB_ROOT  = baSrcPaper/tbb
SUBDIRS   = randmat thresh winnow outer product

# Color Codes
BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'


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
	@echo -e making ... $(BLUE)DASH GNU$(NC);\
	if [[ $$CC == *icc ]] ; then module swap intel gnu; fi;\
	$(MAKE) --no-print-directory dash

dash-intel:
	@echo -e making ... $(BLUE)DASH Intel$(NC);\
	if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) --no-print-directory dash

dash: MAKECMDGOALS = all
dash: $(SUBDIRS)

dash-clean: MAKECMDGOALS = clean
dash-clean: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --no-print-directory -C $@ $(MAKECMDGOALS)
####################    </DASH>    ##############################
####################     <GO>      ##############################
.PHONY: go go-clean $(SUBDIRS:%=go-%) $(SUBDIRS:%=go-clean-%)
#### ------------------------------------------------------------
go-clean: $(SUBDIRS:%=go-clean-%) go-clean-chain
$(SUBDIRS:%=go-clean-%) go-clean-chain:
	-rm $(GO_ROOT)/$(@:go-clean-%=%)/expertpar/main

go:
	@echo -e making ... $(BLUE)GO$(NC)
	@$(MAKE) --no-print-directory go_

go_:$(SUBDIRS:%=go-%) go-chain

$(SUBDIRS:%=go-%) go-chain:
	@$(MAKE) --no-print-directory -C $(GO_ROOT)/$(@:go-%=%)/expertpar main
####################     </GO>     ##############################
####################    <Chapel>   ##############################
.PHONY: chpl chpl_ chpl-clean $(SUBDIRS:%=chpl-%) $(SUBDIRS:%=chpl-clean-%)
#### ------------------------------------------------------------
chpl:
	@echo -e making ... $(BLUE)Chapel$(NC);\
	if [[ $$CC == *icc ]] ; then module swap intel gnu; fi;\
	$(MAKE) --no-print-directory chpl_

chpl_: $(SUBDIRS:%=chpl-%) chpl-chain

$(SUBDIRS:%=chpl-%) chpl-chain:
	@$(MAKE) --no-print-directory -C $(CHPL_ROOT)/$(@:chpl-%=%)/expertpar main

chpl-clean: $(SUBDIRS:%=chpl-clean-%) chpl-clean-chain
$(SUBDIRS:%=chpl-clean-%) chpl-clean-chain:
	-rm $(CHPL_ROOT)/$(@:chpl-clean-%=%)/expertpar/main
####################    </Chapel>  ##############################
####################     <Cilk>    ##############################
.PHONY: cilk cilk_ cilk-clean $(SUBDIRS:%=cilk-%) $(SUBDIRS:%=cilk-clean-%)
#### ------------------------------------------------------------
cilk:
	@echo -e making ... $(BLUE)Cilk$(NC);\
	if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) --no-print-directory cilk_

cilk_: $(SUBDIRS:%=cilk-%) cilk-chain

$(SUBDIRS:%=cilk-%) cilk-chain:
	@$(MAKE) --no-print-directory -C $(CILK_ROOT)/$(@:cilk-%=%)/expertpar main

cilk-clean: $(SUBDIRS:%=cilk-clean-%) cilk-clean-chain

cilk-clean-chain:
	-rm $(CILK_ROOT)/chain/expertpar/main
	-rm $(CILK_ROOT)/chain/expertpar/randmat.o
	-rm $(CILK_ROOT)/chain/expertpar/thresh.o
	-rm $(CILK_ROOT)/chain/expertpar/winnow.o
	-rm $(CILK_ROOT)/chain/expertpar/outer.o
	-rm $(CILK_ROOT)/chain/expertpar/product.o

$(SUBDIRS:%=cilk-clean-%):
	-rm $(CILK_ROOT)/$(@:cilk-clean-%=%)/expertpar/main
####################    </Cilk>    ##############################
####################     <TBB>     ##############################
.PHONY: tbb tbb_ tbb-clean $(SUBDIRS:%=tbb-%) $(SUBDIRS:%=tbb-clean-%)
#### ------------------------------------------------------------
tbb:
	@echo -e making ... $(BLUE)TBB$(NC);\
	if [[ $$CC == *gcc ]] ; then module swap gnu intel; fi;\
	$(MAKE) --no-print-directory tbb_

tbb_: $(SUBDIRS:%=tbb-%) tbb-chain

$(SUBDIRS:%=tbb-%) tbb-chain:
	@$(MAKE) --no-print-directory -C $(TBB_ROOT)/$(@:tbb-%=%)/expertpar main

tbb-clean: $(SUBDIRS:%=tbb-clean-%) tbb-clean-chain

tbb-clean-chain:
	-rm $(TBB_ROOT)/chain/expertpar/main
	-rm $(TBB_ROOT)/chain/expertpar/randmat.o
	-rm $(TBB_ROOT)/chain/expertpar/thresh.o
	-rm $(TBB_ROOT)/chain/expertpar/winnow.o
	-rm $(TBB_ROOT)/chain/expertpar/outer.o
	-rm $(TBB_ROOT)/chain/expertpar/product.o

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
	echo "DASH ..... run with: nRowsCols:$(nRowsCols), seed:$(seed), thresh:$(thresh), winnow_nelts:$(winnow_nelts)";\
	$(MAKE) --no-print-directory dash-run_ ;\
	echo "DASH ..... finished"

dash-run_:
	@echo $(nRowsCols) $(nRowsCols) > run_results/dash_randmat_out
	@echo $(nRowsCols) $(nRowsCols) $(seed) | ./randmat/randmat >> run_results/dash_randmat_out
	@{ cat run_results/dash_randmat_out; echo $(thresh); } | ./thresh/thresh > run_results/dash_thresh_out
#@x=$$(<run_results/dash_thresh_out); y="$${x//[^1]}"; echo $${#y} >> run_results/dash_thresh_out;	#count the '1's in the thresh output and add the number to the output
	@echo $(winnow_nelts) >> run_results/dash_thresh_out
	@cat run_results/dash_randmat_out run_results/dash_thresh_out | ./winnow/main > run_results/dash_winnow_out
	@./outer/outer < run_results/dash_winnow_out > run_results/dash_outer_out
# @echo $(thresh) >> run_results/dash_randmat_out
# @./thresh/thresh < run_results/dash_randmat_out > run_results/dash_thresh_out


.PHONY: chpl-run chpl-run_ chpl-run-randmat
chpl-run:| run_results
	@if [[ $$CC == *icc ]] ; then module swap intel gnu &> /dev/null; fi;\
	echo "Chapel ... run with: nRowsCols:$(nRowsCols), seed:$(seed), thresh:$(thresh), winnow_nelts:$(winnow_nelts)";\
	$(MAKE) --no-print-directory chpl-run_ ;\
	echo "Chapel ... finished"

chpl-run_:
	@echo $(nRowsCols) $(nRowsCols) > run_results/chpl_randmat_out
	@echo $(nRowsCols) $(nRowsCols) $(seed) | $(CHPL_ROOT)/randmat/expertpar/main >> run_results/chpl_randmat_out
	@{ cat run_results/chpl_randmat_out; echo $(thresh); } | $(CHPL_ROOT)/thresh/expertpar/main > run_results/chpl_thresh_out
#@x=$$(<run_results/chpl_thresh_out); y="$${x//[^1]}"; echo $${#y} >> run_results/chpl_thresh_out;	#count the '1's in the thresh output and add the number to the output
	@echo $(winnow_nelts) >> run_results/chpl_thresh_out
	@cat run_results/chpl_randmat_out run_results/chpl_thresh_out | $(CHPL_ROOT)/winnow/expertpar/main > run_results/chpl_winnow_out
	@$(CHPL_ROOT)/outer/expertpar/main < run_results/chpl_winnow_out > run_results/chpl_outer_out
	@$(CHPL_ROOT)/product/expertpar/main < run_results/chpl_outer_out > run_results/chpl_product_out
	@echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(CHPL_ROOT)/chain/expertpar/main > run_results/chpl_chain_out


.PHONY: go-run go-run_ go-run-randmat
go-run:| run_results
	@echo "Go ....... run with: nRowsCols:$(nRowsCols), seed:$(seed), thresh:$(thresh), winnow_nelts:$(winnow_nelts)"
	@$(MAKE) --no-print-directory go-run_
	@echo "Go ....... finished"

go-run_:
	@echo $(nRowsCols) $(nRowsCols) > run_results/go_randmat_out
	@echo $(nRowsCols) $(nRowsCols) $(seed) | $(GO_ROOT)/randmat/expertpar/main >> run_results/go_randmat_out
	@{ cat run_results/go_randmat_out; echo $(thresh); } | $(GO_ROOT)/thresh/expertpar/main > run_results/go_thresh_out
#@x=$$(<run_results/go_thresh_out); y="$${x//[^1]}"; echo $${#y} >> run_results/go_thresh_out;	#count the '1's in the thresh output and add the number to the output
	@echo $(winnow_nelts) >> run_results/go_thresh_out
	@cat run_results/go_randmat_out run_results/go_thresh_out | $(GO_ROOT)/winnow/expertpar/main > run_results/go_winnow_out
	@$(GO_ROOT)/outer/expertpar/main < run_results/go_winnow_out > run_results/go_outer_out
	@$(GO_ROOT)/product/expertpar/main < run_results/go_outer_out > run_results/go_product_out
	@echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(GO_ROOT)/chain/expertpar/main > run_results/go_chain_out


.PHONY: cilk-run cilk-run_ cilk-run-randmat
cilk-run: |run_results
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;\
	echo "Cilk ..... run with: nRowsCols:$(nRowsCols), seed:$(seed), thresh:$(thresh), winnow_nelts:$(winnow_nelts)";\
	$(MAKE) --no-print-directory cilk-run_ ;\
	echo "Cilk ..... finished"

cilk-run_:
	@echo $(nRowsCols) $(nRowsCols) > run_results/cilk_randmat_out
	@echo $(nRowsCols) $(nRowsCols) $(seed) | $(CILK_ROOT)/randmat/expertpar/main >> run_results/cilk_randmat_out
	@{ cat run_results/cilk_randmat_out; echo $(thresh); } | $(CILK_ROOT)/thresh/expertpar/main > run_results/cilk_thresh_out
#@x=$$(<run_results/cilk_thresh_out); y="$${x//[^1]}"; echo $${#y} >> run_results/cilk_thresh_out;	#count the '1's in the thresh output and add the number to the output
	@echo $(winnow_nelts) >> run_results/cilk_thresh_out
	@cat run_results/cilk_randmat_out run_results/cilk_thresh_out | $(CILK_ROOT)/winnow/expertpar/main > run_results/cilk_winnow_out
	@$(CILK_ROOT)/outer/expertpar/main < run_results/cilk_winnow_out > run_results/cilk_outer_out
	@$(CILK_ROOT)/product/expertpar/main < run_results/cilk_outer_out > run_results/cilk_product_out
	@echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(CILK_ROOT)/chain/expertpar/main > run_results/cilk_chain_out


.PHONY: tbb-run tbb-run_ tbb-run-randmat
tbb-run:| run_results
	@if [[ $$CC == *gcc ]] ; then module swap gnu intel &> /dev/null; fi;\
	echo "TBB ...... run with: nRowsCols:$(nRowsCols), seed:$(seed), thresh:$(thresh), winnow_nelts:$(winnow_nelts)";\
	$(MAKE) --no-print-directory tbb-run_ ;\
	echo "TBB ...... finished"

tbb-run_:
	@echo $(nRowsCols) $(nRowsCols) > run_results/tbb_randmat_out
	@echo $(nRowsCols) $(nRowsCols) $(seed) | $(TBB_ROOT)/randmat/expertpar/main >> run_results/tbb_randmat_out
	@{ cat run_results/tbb_randmat_out; echo $(thresh); } | $(TBB_ROOT)/thresh/expertpar/main > run_results/tbb_thresh_out
#@x=$$(<run_results/tbb_thresh_out); y="$${x//[^1]}"; echo $${#y} >> run_results/tbb_thresh_out;	#count the '1's in the thresh output and add the number to the output
	@echo $(winnow_nelts) >> run_results/tbb_thresh_out
	@cat run_results/tbb_randmat_out run_results/tbb_thresh_out | $(TBB_ROOT)/winnow/expertpar/main > run_results/tbb_winnow_out
	@$(TBB_ROOT)/outer/expertpar/main < run_results/tbb_winnow_out > run_results/tbb_outer_out
	@$(TBB_ROOT)/product/expertpar/main < run_results/tbb_outer_out > run_results/tbb_product_out
	@echo $(nRowsCols) $(seed) $(thresh) $(winnow_nelts) | $(TBB_ROOT)/chain/expertpar/main > run_results/tbb_chain_out



SUCCESS:=$(GREEN)SUCCESS$(NC)
FAILURE:=$(RED)FAILURE$(NC)


#export nRowsCols=8 seed=$$(( 1 + RANDOM % 666 )) thresh=33 winnow_nelts=10;\
#export nRowsCols=$$vNRowsCols seed=$$(( 1 + RANDOM % 666 )) thresh=$$vThresh winnow_nelts=$$vWinnow_nelts;\

.PHONY: test
test:
	@echo Enter start nRowsCols and then max nRowsCols:; read start; read i;\
	echo ---------- ;\
	for (( succ=0 , fail=0 , subsucc=0, warnings=0, globWarnings=0, vNRowsCols=$$start; vNRowsCols <= $$i ; ++vNRowsCols )) ; do \
	for (( vThresh=0; vThresh<= 100; vThresh+=1 )) ;                   do \
	for (( vWinnow_nelts=2; vWinnow_nelts<= 15 ; vWinnow_nelts+=1 )) ; do \
	export nRowsCols=15 seed=197 thresh=36 winnow_nelts=10;\
		$(MAKE) --no-print-directory run-all;\
		echo comparing ... ;\
		subsucc=0;\
		if [[ $$(diff -swB run_results/tbb_randmat_out run_results/chpl_randmat_out) == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Chapel randmat output differs from TBB!"$(NC); fi;\
		if [[ $$(diff -swB run_results/tbb_randmat_out run_results/go_randmat_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"GO randmat output differs from TBB!"$(NC);     fi;\
		if [[ $$(diff -swB run_results/tbb_randmat_out run_results/cilk_randmat_out) == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Cilk randmat output differs from TBB!"$(NC);   fi;\
		if [[ $$(diff -swB run_results/tbb_randmat_out run_results/dash_randmat_out) == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"DASH randmat output differs from TBB!"$(NC);   fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)"randmat comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at randmat comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		subsucc=0;\
		if [[ $$(diff -swB run_results/tbb_thresh_out run_results/go_thresh_out)     == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"GO thresh output differs from TBB!"$(NC);     fi;\
		if [[ $$(diff -swB run_results/tbb_thresh_out run_results/cilk_thresh_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Cilk thresh output differs from TBB!"$(NC);   fi;\
		if [[ $$(diff -swB run_results/tbb_thresh_out run_results/chpl_thresh_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Chapel thresh output differs from TBB!"$(NC); fi;\
		if [[ $$(diff -swB run_results/tbb_thresh_out run_results/dash_thresh_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"DASH thresh output differs from TBB!"$(NC);   fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)"thresh  comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at thresh comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		subsucc=0;\
		if [[ $$(diff -swB run_results/tbb_winnow_out run_results/go_winnow_out)     == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"GO winnow output differs from TBB!"$(NC);     fi;\
		if [[ $$(diff -swB run_results/tbb_winnow_out run_results/cilk_winnow_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Cilk winnow output differs from TBB!"$(NC);   fi;\
		if [[ $$(diff -swB run_results/tbb_winnow_out run_results/chpl_winnow_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"Chapel winnow output differs from TBB!"$(NC); fi;\
		if [[ $$(diff -swB run_results/tbb_winnow_out run_results/dash_winnow_out)   == *are\ identical* ]] ;then (( ++subsucc )); else echo -e $(RED)"DASH winnow output differs from TBB!"$(NC);   fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)"winnow  comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at winnow comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		subsucc=0;\
		if [[     $$(        diff -swB run_results/tbb_outer_out run_results/go_outer_out)       == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_outer_out run_results/go_outer_out)       == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_outer_out run_results/go_outer_out)       == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"GO outer output differs from TBB!"$(NC);     fi;\
		if [[     $$(        diff -swB run_results/tbb_outer_out run_results/cilk_outer_out)     == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_outer_out run_results/cilk_outer_out)     == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_outer_out run_results/cilk_outer_out)        == *are\ equal*  ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Cilk outer output differs from TBB!"$(NC);   fi;\
		if [[     $$(        diff -swB run_results/tbb_outer_out run_results/chpl_outer_out)     == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_outer_out run_results/chpl_outer_out)     == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_outer_out run_results/chpl_outer_out)     == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Chapel outer output differs from TBB!"$(NC); fi;\
		if [[     $$(        diff -swB run_results/tbb_outer_out run_results/dash_outer_out)     == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_outer_out run_results/dash_outer_out)     == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_outer_out run_results/dash_outer_out)     == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"DASH outer output differs from TBB!"$(NC); fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)"outer   comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at outer comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		subsucc=0;\
		if [[     $$(        diff -swB run_results/tbb_product_out run_results/go_product_out)   == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_product_out run_results/go_product_out)   == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_product_out run_results/go_product_out)   == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"GO product output differs from TBB!"$(NC);     fi;\
		if [[     $$(        diff -swB run_results/tbb_product_out run_results/cilk_product_out) == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_product_out run_results/cilk_product_out) == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_product_out run_results/cilk_product_out) == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Cilk product output differs from TBB!"$(NC);   fi;\
		if [[     $$(        diff -swB run_results/tbb_product_out run_results/chpl_product_out) == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_product_out run_results/chpl_product_out) == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_product_out run_results/chpl_product_out) == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Chapel product output differs from TBB!"$(NC); fi;\
		case $$subsucc in \
			3) (( ++succ )) ; echo -e $(GREEN)"product comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at product comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		subsucc=0;\
		if [[     $$(        diff -swB run_results/tbb_product_out run_results/tbb_chain_out)   == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_product_out run_results/tbb_chain_out)   == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_product_out run_results/tbb_chain_out)   == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"TBB chain output differs from TBB Product Output!"$(NC);  fi;\
		if [[     $$(        diff -swB run_results/tbb_chain_out run_results/go_chain_out)      == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_chain_out run_results/go_chain_out)      == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_chain_out run_results/go_chain_out)      == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"GO chain output differs from TBB chain Output!"$(NC);     fi;\
		if [[     $$(        diff -swB run_results/tbb_chain_out run_results/cilk_chain_out)    == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_chain_out run_results/cilk_chain_out)    == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_chain_out run_results/cilk_chain_out)    == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Cilk chain output differs from TBB chain Output!"$(NC);   fi;\
		if [[     $$(        diff -swB run_results/tbb_chain_out run_results/chpl_chain_out)    == *are\ identical* ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.9e-2 run_results/tbb_chain_out run_results/chpl_chain_out)    == *are\ equal*     ]] ;then (( ++subsucc ));\
			elif [[ $$(numdiff -a 0.1    run_results/tbb_chain_out run_results/chpl_chain_out)    == *are\ equal*     ]] ;then (( ++subsucc )); (( ++warnings ));\
			else echo -e $(RED)"Chapel chain output differs from TBB chain Output!"$(NC); fi;\
		case $$subsucc in \
			4) (( ++succ )) ; echo -e $(GREEN)"chain   comparison without differences."$(NC);;\
			*) (( ++fail )) ; echo -e $(RED)"error at chain comparison with: "$(NC) "nRowsCols:" $$nRowsCols " seed:" $$seed " thresh:" $$thresh " winnow_nelts:" $$winnow_nelts >> error.log ;;\
		esac;\
		if (( warnings > 0 )) ;then\
			echo -e $(BLUE)warnings at: $(NC) nRowsCols: $$nRowsCols  seed: $$seed thresh: $$thresh winnow_nelts: $$winnow_nelts >> error.log ;\
			(( globWarnings += warnings )); warnings=0; fi;\
		echo --------------------------------------;\
		break;\
	done ;\
	break;\
	 done ;\
	 break;\
	  done ;\
	echo -e $(GREEN)successes$(NC):$$succ $(RED)failures$(NC):$$fail $(BLUE)warnings$(NC):$$globWarnings;\

.PHONY: run-clean
run-clean:
	-rm -r run_results

################  </Validation Test>   ##########################
