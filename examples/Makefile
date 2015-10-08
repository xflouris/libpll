BINDIR=bin
SUBDIRS=heterotachy \
	mixture \
	model-optimize \
	newick-fasta-rooted \
	newick-fasta-unrooted \
	partial-traversal \
	protein-list \
	rooted \
	rooted-tacg \
	unrooted

all:
	@mkdir -p $(BINDIR)
	@for a in $(SUBDIRS); do \
	    if [ -d $$a ]; then \
            echo "processing test $$a"; \
            $(MAKE) -C $$a; \
	    cp $$a/$$a bin/; \
        fi; \
	done;
	@echo "Done!"

subdir:
	@mkdir -p $(BINDIR)
	@if [ -d $(subdir) ]; then \
            echo "processing test $(subdir)"; \
            $(MAKE) -C $(subdir); \
	    cp $(subdir)/$(subdir) bin/; \
	else \
	    echo "ignoring $(subdir); \
        fi;

clean:
	@for a in $(SUBDIRS); do \
	    if [ -d $$a ]; then \
            echo "processing test $$a"; \
            $(MAKE) -C $$a clean; \
        fi; \
	done;
	@rm -rf $(BINDIR)
