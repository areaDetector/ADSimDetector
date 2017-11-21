#FILENAME:      Makefile
#DESCRIPTION:   convert MEDM screens to other GUI formats


SUBDIRS = ui

.PHONY: all subdirs $(SUBDIRS)

all :: subdirs

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
