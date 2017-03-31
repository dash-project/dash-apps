include ./make.defs

SUBDIRS = randmat thresh winnow outer product

.PHONY: all $(SUBDIRS)
all: $(SUBDIRS)

clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

