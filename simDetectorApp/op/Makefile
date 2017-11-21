#FILENAME:      Makefile
#DESCRIPTION:   convert MEDM screens to other GUI formats


# add additional directories once the Makefiles have been created
# TODO: switch to system of EPICS Makefile rules
SUBDIRS = ui edl


.PHONY: all subdirs $(SUBDIRS)

all :: subdirs

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
