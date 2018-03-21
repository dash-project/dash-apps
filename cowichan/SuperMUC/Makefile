include ./make.defs

SUBDIRS   = randmat thresh winnow outer product chain


.PHONY: all clean all-clean
#### ------------------------------------------------------------
all: dash

clean: dash-clean
all-clean: clean
####################  </Top Level Targets>  #####################
####################       <DASH>           #####################
.PHONY: dash dash-clean $(SUBDIRS)
#### ------------------------------------------------------------

dash: MAKECMDGOALS = all
dash: $(SUBDIRS)

dash-clean: MAKECMDGOALS = clean
dash-clean: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) --no-print-directory -C $@ $(MAKECMDGOALS)
####################    </DASH>    ##############################

