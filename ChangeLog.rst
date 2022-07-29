=======
Changes
=======
todo in future versions:
- Supporting for *de novo* GTF file (no annotated start and stop codons).

v1.2.15 (2022.7.29)
-----------------
- Fixed a small bug caused by None values in levels ("metaplots.py", level = list(sorted(levels))[0] line 49)

v1.2.14 (2022.1.25)
-----------------
- Add "--dependence_test" option for measuring the dependence of densities between frame 1 and frame 2 of each ORF. This test would help determine whether the adjustment should be performed on the combined p-values.

v1.2.13 (2021.12.26)
-----------------
- Add "--stouffer_adj" option for adjusting cominbed p-values to account for the dependence between two tests (i.e. F0 vs F1 and F0 vs F2)
- Add "--pval_adj" option to correct p-values for multiple testing
- Report alternative annotations of each ORF based on the transcripts other than the longest one (some ORFs in various transcript isoforms but from same gene might be annotated as different types)

v1.2.12 (2021.10.28)
-----------------
- Fixed a small bug caused by an update of h5py; Added a new parameter "--plot-annotated-orf".

v1.2.11 (2018.12.20)
-----------------
- Several small updates.

v1.2.10 (2018.3.25)
-----------------
- Add "-i" option to metaplots command for multiple input files
- Supporting for non-stranded library sequencing
- Update the document

v1.2.9 (2018.3.5)
-----------------
- Small bug fixed in metaplots.
- Added more test function.
- Change the default backend of matplotlib.

v1.2.8 (2017.12.30)
-------------------
- Supporting for Python3
- Small changes in test_func.py
- Update metaplot script

v1.2.7 (2017.7.25)
------------------
- Add "--version" option to show the version number of package.
- Add RiboCode_onestep command which performs entire steps for detecting ORFs.
- Supporting for custom GTF file

v1.2.6 (2017.5.17)
------------------
- Fixed a small bug in ORFcount command.
- Included package into bioconda.

v1.2.5 (2017.5.12)
------------------
- Add support for counting RPF reads aligned to predicted ORFs.
- Update document.

v1.2.4 (2017.5.4)
-----------------
- Add support for outputting ORF results in gtf format.
- Fix a bug where some ORFs' genomic coordinates are wrong.
- Other optimizations on code and documents.

v1.2.3 (2017.4.12)
------------------
- Automatically selecting periodic ribo-seq read lengths and determining the P-site location for each read length.
- Fixed a small bug.

v1.2.0 (2017.3.8)
-----------------
- version 1.2.0 release !
