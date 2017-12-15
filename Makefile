

.PHONY : run nbhtml x3d

nbhtml : $(shell python get_notebooknames.py "../temp/nb_html/" ".html")

arcusversionfile = $(shell python -c "from __future__ import print_function; import arcus.version as vers; print(vers.__file__)")

.PRECIOUS : %.processed.ipynb

x3d : plot_scripts/plot_master.py $(arcusversionfile)
	mkdir -p ../temp/x3d
	python plot_scripts/plot_master.py ../temp/x3d

../temp/nb_processed/%.processed.ipynb : notebooks/%.ipynb $(arcusversionfile)
	mkdir -p ../temp/nb_processed
	jupyter nbconvert --ExecutePreprocessor.timeout=1800 --to notebook --execute $< --output ../$@

../temp/nb_html/%.html : ../temp/nb_processed/%.processed.ipynb
	mkdir -p ../temp/nb_html
	jupyter nbconvert --to html $< --output ../$@

website: nbhtml x3d
	mkdir -p ../web_out
	python website/build_site.py
