

.PHONY : run nbhtml nbpdf x3d

nbhtml : $(shell python get_notebooknames.py "nb_proc/" ".html")
nbpdf : $(shell python get_notebooknames.py "nb_proc/" ".pdf")

arcusversionfile = $(shell python -c "import arcus.version as vers; print(vers.__file__)")

.PRECIOUS : nb_proc/%.processed.ipynb

x3d : plot_scripts/plot_master.py $(arcusversionfile)
	mkdir -p ../temp/x3d
	python plot_scripts/plot_master.py ../temp/x3d

nb_proc/%.processed.ipynb : notebooks/%.ipynb $(arcusversionfile)
	mkdir -p nb_proc
	# Set kernel name to how the kernel is called when inside the environment
	# Notebook meta data may have global kernel name
	# Note output dir is relative to input dir
	# so we build in same dir and use ../$@
	cd notebooks && jupyter nbconvert --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=1800 --to notebook --execute $(notdir $<) --output ../$@

nb_proc/%.html : nb_proc/%.processed.ipynb
	cd nb_proc && jupyter nbconvert --to html $(notdir $<) --output $(notdir $@)

nb_proc/%.pdf : nb_proc/%.processed.ipynb
	cd nb_proc && jupyter nbconvert --no-input --to PDF $(notdir $<) --output $(notdir $@)

website: nbhtml x3d nbpdf
	mkdir -p ../web_out
	python website/build_site.py
