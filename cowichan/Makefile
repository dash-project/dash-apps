include ./make.defs

SUBDIRS = randmat thresh winnow outer product

.PHONY: all clean $(SUBDIRS)
all: $(SUBDIRS)

clean: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

