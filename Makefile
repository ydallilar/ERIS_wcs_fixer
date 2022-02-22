
DESTDIR=/usr/local/bin

install:
	install eris_wcs_fixer.py $(DESTDIR)/

clean:
	rm -rf $(DESTDIR)/eris_wcs_fixer.py

