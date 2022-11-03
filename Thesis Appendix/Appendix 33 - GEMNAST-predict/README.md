# GEMNAST-predict

The current prototype utilises abundance data from a given patient (specified by the user) and determines the percentage (relative abundance) of primary degraders (strains that can degrade such fibre and release monosaccharides) and secondary degraders (strains that can utilise released monosaccharides to export a especific intermediate metabolite)

The user has to provide three main inputs:
- Patient ID
- Fibre (select one from corresponding file from Data_files folder)
- Intermediate metabolite (refer to script for options)

The script utilises GEMNAST-generated data to identify primary and secondary degraders. The required files are provided in the Data_files folder.

Additionally, patient abundance data is also provided in the Data_files folder.

GEMNAST-predict was written in the Python coding language and requires installation of the latest version of Python. A Python-compatible IDE is also recommended (e.g. PyCharm).

Please refer to in-script commentaries for additional details.
