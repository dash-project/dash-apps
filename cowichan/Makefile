include ./make.defs

SUBDIRS = randmat thresh winnow outer product

.PHONY: all $(SUBDIRS)
#all: subdirs
all: $(SUBDIRS)

#.PHONY: subdirs $(SUBDIRS)

#subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

#.PHONY: clean
#clean: $(SUBDIRS)
#	$(MAKE) -C $? clean
