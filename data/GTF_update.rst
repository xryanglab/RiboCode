=============================
Preparing the GTF format file
=============================

The current RiboCode version supports the standard GTF format, and it also requires
the GTF file to be satisfied with three-level hierarchy annotations (genes contain transcripts that contain exons
and optionally, a CDS). This type of file can be obtained from ENSEMBL/GENCODE database. Those from other source or
the custom GTF file may lack the gene and transcript annotation information. RiboCode provide a command "GTFupdate" which
can add these two information for GTF file:

.. code-block:: bash

  GTFupdate original.gtf > updated.gtf

Example:
--------

**original GTF:**

1 unknown exon 11874 12227 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";

1 unknown exon 12613 12721 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";

1 unknown exon 13221 14409 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";

**after updating:**

*1 unknown gene 11874 14409 . + . gene_id "DDX11L1"; gene_name "DDX11L1";*

*1 unknown transcript 11874 14409 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";*

1 unknown exon 11874 12227 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";

1 unknown exon 12613 12721 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";

1 unknown exon 13221 14409 . + . gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018.2";


the standard GTF format
-----------------------
see detail in this website: https://en.wikipedia.org/wiki/GENCODE

+---------------+-----------------------+
| column-number +         content       |
+===============+=======================+
|      1        +       chromosome      |
+---------------+-----------------------+
|      2        +    source (not used)  |
+---------------+-----------------------+
|      3        +  feature type         |
+---------------+-----------------------+
|      4        + genomic start location|
+---------------+-----------------------+
|      5        + genomic end location  |
+---------------+-----------------------+
|      6        +    score(not used)    |
+---------------+-----------------------+
|      7        +   genomic strand      |
+---------------+-----------------------+
|      8        +genomic phase(not used)|
+---------------+-----------------------+
|      9        + attributes (see below)|
+---------------+-----------------------+

Description of attributes in 9th column of the GTF file
-------------------------------------------------------
- gene_id
- transcript_id
- gene_type (optional)
- gene_name  (optional)
- transcript_id  (optional)
- transcript_type  (optional)
- level (optional)
- cdds  (optional)

