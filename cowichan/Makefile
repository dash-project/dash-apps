include ./make.defs

SUBDIRS = randmat thresh winnow outer product
GOROOT = baSrcPaper/go

PAPERSUB = expertpar

.PHONY: all clean $(SUBDIRS) printenv all-clean

all:   $(SUBDIRS)
clean: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

#postfix $(SUBDIRS:%=%-go)
.PHONY: go go-clean $(SUBDIRS:%=go-%) $(SUBDIRS:%=go-clean-%)
#go:MAKECMDGOALS=main
go:$(SUBDIRS:%=go-%)

go-clean: $(SUBDIRS:%=go-clean-%)

$(SUBDIRS:%=go-clean-%):
	rm $(GOROOT)/$(@:go-clean-%=%)/expertpar/main

$(SUBDIRS:%=go-%): #setpapersub
	$(MAKE) -C $(GOROOT)/$(@:go-%=%)/$(PAPERSUB) main

#setpapersub:
#	PAPERSUB=$(MAKECMDGOALS)
