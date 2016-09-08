dist_man_MANS=libpll.3

if HAVE_PS2PDF

doc_DATA = libpll_manual.pdf

libpll_manual.pdf : libpll.3
	TEMP=$$(mktemp temp.XXXXXXXX) ; \
	if [ $$(uname) == "Darwin" ] ; then \
		${SED} -e 's/\\-/-/g' $< | \
		iconv -f UTF-8 -t ISO-8859-1 > $$TEMP ; \
	else \
		${SED} -e 's/\\-/-/g' $< > $$TEMP ; \
	fi ; \
	man -t ./$$TEMP | ps2pdf -sPAPERSIZE=a4 - $@ ; \
	rm $$TEMP

CLEANFILES=libpll_manual.pdf

endif
