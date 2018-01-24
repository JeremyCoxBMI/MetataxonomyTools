# MetataxonomyTools
Libraries to help process Metataxonomy data, specifically NCBI gid, accession_numbers, and taxon ID.

I have taken all the scripts I have written and created libraries for my own use.  And why not share?


## Features
  * I use these libraries to parse alignments against metagenome databases, for example results from IMSA-A (my project on github).  
  ### Alignment files
  * Accepts alignment files blast format6 (BLAST) or sam (BWA, BOWTIE2, etc) or psl (BLAT).
  * Generates useful alignment statistics for all file types; enables comparative analyses.
  * Use Genbank IDs (gid) or NCBI Accession numbers.  Note that NCBI has deprecated the use of gid.
  ### Taxonomy
  * Uses local data files improving speed over Eutils. 
  * Convert taxon ID to other taxons in the lineage
  * Create cladograms
  
## Dependencies
  * Written for Python 2.7.  
  * I recommend using Anaconda when running on Windows.
  * ETE2 (etetoolkit.org) :: for making cladograms only

## Versions
  ### v0.0  
  Project creation 2018-01-24
  README.md drafted
