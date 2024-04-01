import bibtexparser
# version < 2.0.0


with open("report.bib") as bibtex_file:
    bib_database = bibtexparser.load(bibtex_file)

for entry in bib_database.entries:
    if "doi" in entry:
        entry["note"] = f'[{entry["doi"]}]'

with open("report_with_doi.bib", "w") as bibtex_file:
    bibtexparser.dump(bib_database, bibtex_file)