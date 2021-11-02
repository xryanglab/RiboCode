====================================================
Detect translated ORFs using ribosome-profiling data
====================================================

|BuildStatus| |PyPI| |PythonVersions| |BioConda| |Publish| |downloads|

*RiboCode* is a very simple but high-quality computational algorithm to
identify genome-wide translated ORFs using ribosome-profiling data.

Dependencies:
-------------

- pysam

- pyfasta

- h5py

- Biopython

- Numpy

- Scipy

- matplotlib

- HTSeq

Installation
------------

*RiboCode* can be installed like any other Python packages. Here are some popular ways:

* Install via pypi:

.. code-block:: bash

  pip install ribocode

* Install via conda:

.. code-block:: bash

  conda install -c bioconda ribocode

* Install from source:

.. code-block:: bash

  git clone https://www.github.com/xzt41/RiboCode
  cd RiboCode
  python setup.py install

* Install from local:

.. code-block:: bash

  pip install RiboCode-*.tar.gz

If you have not administrator permission, you need to install *RiboCode* locally in you own directory by adding the
option ``--user`` in the above command. Then, you need to define ``~/.local/bin/`` in ``PATH`` variable,
and ``~/.local/lib/`` in ``PYTHONPATH`` variable. For example, if you are using the bash shell, you should add the following lines to your ``~/.bashrc`` file:

.. code-block:: bash

  export PATH=$PATH:$HOME/.local/bin/
  export PYTHONPATH=$HOME/.local/lib/python2.7

then, source your ``~/.bashrc`` file using this command:

.. code-block:: bash

  source ~/.bashrc

Users can also update or uninstall package through one of the following commands:

.. code-block:: bash

  pip install --upgrade RiboCode # upgrade
  pip uninstall RiboCode # uninstall
  conda update -c bioconda ribocode # upgrade
  conda remove ribocode # uninstall

Tutorial to analyze ribosome-profiling data and run *RiboCode*
--------------------------------------------------------------

Here, we use the `HEK293 dataset`_ as an example to illustrate the use of *RiboCode* and demonstrate typical workflow.
Please make sure the path and file name are correct.

1. **Required files** 

   The genome FASTA file, GTF file for annotation can be downloaded from:


   http://www.gencodegenes.org

   or from:

   http://asia.ensembl.org/info/data/ftp/index.html

   http://useast.ensembl.org/info/data/ftp/index.html


   For example, the required files in this tutorial can be downloaded from following URL:

   GTF: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

   FASTA: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

   |Important| The GTF file required by *RiboCode* should include three-level hierarchy
   annotations: genes,transcripts and exons. Some GTF files may lack the gene and transcript
   annotations, users can added these annotations using the "GTFupdate" command in *RiboCode*.
   Please refer to `GTF_update.rst`_ for more information.

   The raw Ribo-seq FASTQ file can be downloaded using fastq-dump tool from `SRA_Toolkit`_:

   .. code-block:: bash

      fastq-dump -A <SRR1630831>

2. **Trimming adapter sequence for ribo-seq data**

   Using cutadapt program https://cutadapt.readthedocs.io/en/stable/installation.html

   Example:

   .. code-block:: bash

      cutadapt -m 20 --match-read-wildcards -a (Adapter sequence) -o <Trimmed fastq file> <Input fastq file>


   Here, the adapter sequences for this data had already been trimmed off, so we can skip this step.

3. **Removing ribosomal RNA(rRNA) derived reads**

   Removing rRNA contamination by aligning the trimmed reads to rRNA sequences using `Bowtie`_,
   then keeping the unaligned reads for the next step.

   rRNA sequences are provided in `rRNA.fa`_ file.

   Example:

   .. code-block:: bash

      bowtie-build <rRNA.fa> rRNA
      bowtie -p 8 -norc --un <un_aligned.fastq> -q <SRR1630831.fastq> rRNA <HEK293_rRNA.align>

4. **Aligning the clean reads to reference genome**

   Using STAR program: https://github.com/alexdobin/STAR

   Example:

   (1). Build index

   .. code-block:: bash

      STAR --runThreadN 8 --runMode genomeGenerate --genomeDir <hg19_STARindex>
      --genomeFastaFiles <hg19_genome.fa> --sjdbGTFfile <gencode.v19.annotation.gtf>

   .. _STAR:

   (2). Alignment:

   .. code-block:: bash

      STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir <hg19_STARindex>
      --readFilesIn <un_aligned.fastq>  --outFileNamePrefix <HEK293> --outSAMtype BAM
      SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1
      --outFilterMatchNmin 16 --alignEndsType EndToEnd

5. **Running RiboCode to identify translated ORFs**

   (1). Preparing the transcripts annotation files:

   .. code-block:: bash

      prepare_transcripts -g <gencode.v19.annotation.gtf> -f <hg19_genome.fa> -o <RiboCode_annot>

   (2). Selecting the length range of the RPF reads and identify the P-site locations:

   .. code-block:: bash

      metaplots -a <RiboCode_annot> -r <HEK293Aligned.toTranscriptome.out.bam>


   This step will generate two files: a PDF file plots the aggregate profiles of the distance from the 5'-end
   of reads to the annotated start codons (or stop codons), which is used for examining the P-site periodicity of RPF reads on CDS regions. The P-site config file, which defines the read lengths with
   strong 3-nt periodicity and the associated P-site locations for each length.  In some cases, user may have multiple bam files to predict ORFs
   together in next step, they can use "-i" argument to specify a text file which contains the names of these bam files (
   one file per line)

   .. _RiboCode:

   (3). Detecting translated ORFs using the ribosome-profiling data:

   .. code-block:: bash

      RiboCode -a <RiboCode_annot> -c <config.txt> -l no -g -o <RiboCode_ORFs_result>


   Using the config file generated by last step to specify the information of the bam file and P-site parameters,
   please refer to the example file `config.txt`_ in data folder. The "gtf" or "bed" format file of predicted ORFs can
   be obtained by adding the "-g" or "-b" argument to this command.

   **Explanation of final result files**

   The *RiboCode* generates two text files:
   The "(output file name).txt" contains the information of all predicted ORFs in each transcript.
   The "(output file name)_collapsed.txt" file combines the ORFs having the same stop codon in different transcript
   isoforms: the one harboring the most upstream in-frame ATG will be kept.

   Some column names of the result file::

    - ORF_ID: The identifier of predicated ORF.
    - ORF_type: The type of predicted ORF, which is annotated according to its location to associated CDS. The following ORF categories are reported:

     "annotated" (overlapping with annotated CDS, have the same stop codon with annotated CDS)

     "uORF" (upstream of annotated CDS, not overlapping with annotated CDS)

     "dORF" (downstream of annotated CDS, not overlapping with annotated CDS)

     "Overlap_uORF" (upstream of annotated CDS and overlapping annotated with CDS)

     "Overlap_dORF" (downstream of annotated CDS and overlapping annotated CDS"

     "Internal" (internal ORF of annotated CDS, but in a different reading frame)

     "novel" (from non-coding genes or non-coding transcripts of the coding genes).

    - ORF_tstart, ORF_tstop: the start and end position of ORF relative to its transcript (1-based coordinate)
    - ORF_gstart, ORF_gstop: the start and end position of ORF in the genome (1-based coordinate)
    - pval_frame0_vs_frame1: significance levels of P-site densities of frame0 greater than of frame1
    - pval_frame0_vs_frame2: significance levels of P-site densities of frame0 greater than of frame2
    - pval_combined: integrated P-value by combining pval_frame0_vs_frame1 and pval_frame0_vs_frame2

   **All above three steps can also be easily run by a single command "RiboCode_onestep":**

   .. code-block:: bash

      RiboCode_onestep -g <gencode.v19.annotation.gtf> -f <hg19_genome.fa> -r <HEK293Aligned.toTranscriptome.out.bam>
                       -l no -o <RiboCode_ORFs_result>

   (4). (optional) Plotting the P-sites densities of predicted ORFs

   Using the "plot_orf_density" command, for example:

   .. code-block:: bash

      plot_orf_density -a <RiboCode_annot> -c <config.txt> -t (transcript_id)
      -s (ORF_gstart) -e (ORF_gstop)

   The generated PDF plots can be edited by Adobe Illustrator.

   (5). (optional) Counting the number of RPF reads aligned to ORFs

   The number of reads aligned on each ORF can be counted by the "ORFcount" command which will call the HTSeq-count program.
   Only the reads of a given length will be counted. For those ORF with length longer than a specified value (set by "-e"),
   the RPF reads located in first few and last few codons can be excluded by adjusting the parameters "-f" and "-l".
   For example, the reads with length between 26-34 nt aligned on predicted ORF can be obtained by using below command:

   .. code-block:: bash

      ORFcount -g <RiboCode_ORFs_result.gtf> -r <ribo-seq genomic mapping file> -f 15 -l 5 -e 100 -m 26 -M 34 -o <ORF.counts>

   The reads aligned to first 15 codons and last 5 codons of ORFs and had the length longer than 100 nt will be excluded.
   The "RiboCode_ORFs_result.gtf" file can be generated by `RiboCode`_ command. The "ribo-seq genomic mapping file" is the
   genome-wide mapping file produced by `STAR`_ mapping.


Recipes (FAQ):
--------------
1. **I have a BAM/SAM file aligned to genome, how do I convert it to transcriptome-based mapping file ?**

   You can use STAR aligner to generate the transcriptome-based alignment file by specifying the "--quantMode TranscriptomeSAM" parameters,
   or use the "sam-xlate" command from `UNC Bioinformatics Utilities`_ .

2. **How to use multiple BAM/SAM files to identify ORFs?**

   You can select the read lengths which show strong 3-nt periodicity and the corresponding P-site locations for each
   BAM/SAM file, then list each file and their information in `config.txt`_ file. *RiboCode* will combine the P-site
   densities at each nucleotides of these BAM/SAM files together to predict ORFs.

3. **Generating figures with matplotlib when DISPLAY variable is undefined or invalid**

   When running the "metaplots" or "plot_orf_density" command,  some users received errors similar to the following:

      ``raise RuntimeError('Invalid DISPLAY variable')``

      ``_tkinter.TclError: no display name and no $DISPLAY environment variable``

   The main problem is that default backend of matplotlib is unavailable. The solution is to modify the backend in matplotlibrc file.
   A very simple solution is to set the MPLBACKEND environment variable, either for your current shell or for a single script:

   .. code-block:: bash

      export MPLBACKEND="module:Agg"

   Giving below are non-interactive backends, capable of writing to a file:

      Agg  PS  PDF  SVG  Cairo  GDK

   See also:

   http://matplotlib.org/faq/usage_faq.html#what-is-a-backend

   http://matplotlib.org/users/customizing.html#the-matplotlibrc-file

   http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined


For any questions, please contact:
----------------------------------
Xuerui Yang (yangxuerui@tsinghua.edu.cn)

.. _SRA_Toolkit: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

.. _HEK293 dataset: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1630831

.. _config.txt: https://github.com/xryanglab/RiboCode/blob/master/data/config.txt

.. _rRNA.fa: https://github.com/xryanglab/RiboCode/blob/master/data/rRNA.fa

.. _GTF_update.rst: https://github.com/xryanglab/RiboCode/blob/master/data/GTF_update.rst

.. _UNC Bioinformatics Utilities: https://github.com/mozack/ubu

.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml

.. |PyPI| image:: https://img.shields.io/pypi/v/RiboCode.svg?style=flat-square
   :target: https://pypi.python.org/pypi/RiboCode

.. |PythonVersions| image:: https://img.shields.io/pypi/pyversions/RiboCode.svg?style=flat-square
   :target: https://pypi.python.org/pypi/RiboCode

.. |BioConda| image:: https://img.shields.io/badge/install-bioconda-blue.svg?style=flat-square
   :target: http://bioconda.github.io/recipes/ribocode/README.html

.. |downloads| image:: https://anaconda.org/bioconda/ribocode/badges/downloads.svg
   :target: https://anaconda.org/bioconda/ribocode

.. |Publish| image:: https://img.shields.io/badge/publish-NAR-blue.svg?style=flat-square
   :target: https://doi.org/10.1093/nar/gky179

.. |BuildStatus| image:: https://circleci.com/gh/xryanglab/RiboCode.svg?style=svg
    :target: https://circleci.com/gh/xryanglab/RiboCode

.. |Important| image:: https://img.shields.io/badge/-Note-orange.svg
    :width: 50
    :target: https://github.com/xryanglab/RiboCode/blob/master/data/GTF_update.rst

