
IPYPROC = ../../temp/procnb
HTML = ../../temp/notebook

.PHONY : run html

run : $(IPYPROC)/CATmerrit.processed.ipynb \
$(IPYPROC)/BoomMisalignment.processed.ipynb \
$(IPYPROC)/3vs4sidedboom.processed.ipynb


html :  $(HTML)/CATmerrit.html \
$(HTML)/3vs4sidedboom.html \
$(HTML)/BoomMisalignment.html


.PRECIOUS : %.processed.ipynb

$(IPYPROC)/%.processed.ipynb : %.ipynb
	-mkdir $(IPYPROC)
	jupyter nbconvert --ExecutePreprocessor.timeout=1800 --to notebook --execute $< --output $@

$(HTML)/%.html : $(IPYPROC)/%.processed.ipynb
	-mkdir $(HTML)
	jupyter nbconvert --to html $< --output $@
