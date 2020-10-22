

.PHONY : run nbhtml nbpdf nbipynb x3d

nbhtml : $(shell python get_notebooknames.py "nb_proc/" ".html")
nbpdf : $(shell python get_notebooknames.py "nb_proc/" ".pdf")
nbipynb : $(shell python get_notebooknames.py "nb_proc/" ".ipynb")


Arcusversionfile = $(shell python -c "import arcus.version as vers; print(vers.__file__)")

x3d : plot_scripts/plot_master.py $(arcusversionfile)
	mkdir -p ../temp/x3d
	python plot_scripts/plot_master.py ../temp/x3d

nb_proc/%.ipynb : notebooks/%.ipynb $(arcusversionfile)
	mkdir -p nb_proc
	# Set kernel name to how the kernel is called when inside the environment
	# Notebook meta data may have global kernel name
	# Note output dir is relative to input dir
	# so we build in same dir and use ../$@
	cd notebooks && jupyter nbconvert --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=1800 --to notebook --execute $(notdir $<) --output ../$@

nb_proc/%.html : nb_proc/%.ipynb
	cd nb_proc && jupyter nbconvert --to html $(notdir $<) --output $(notdir $@)

nb_proc/%.pdf : nb_proc/%.ipynb
	cd nb_proc && jupyter nbconvert --no-input --to PDF $(notdir $<) --output $(notdir $@)
