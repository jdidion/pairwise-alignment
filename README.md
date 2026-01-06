# Pairwise alignment

## Algorithms

### General

* Alignment-free tools for pairwise sequence comparison:
    * http://www.combio.pl/alfree/tools/
    * http://www.combio.pl/alfree
    * http://bioinformatics.org.au/tools/decaf+py/
    * Idea: Train a NN to approximate the pairwise alignment distance https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty887/5140215?redirectedFrom=fulltext
    * Optimized embedding function: https://www.biorxiv.org/content/10.1101/2020.05.24.113852v2
* EDLIB: edit distance-based alignment: https://github.com/Martinsos/edlib
* The SeqAn C++ library implements several alignment algorithms:
http://www.sciencedirect.com/science/article/pii/S0168165617315420
* Biotite implementation of optimal alignment: https://www.biotite-python.org/apidoc/biotite.sequence.align.align_optimal.html#biotite.sequence.align.align_optimal
* GASAL2: GPU-accelerated alignment https://github.com/nahmedraja/GASAL2
* Stdaln: Fast dynamic programming algorithm by Heng Li: https://github.com/jts/nanopolish/tree/master/src/thirdparty
* Heng Li review of dynamic programming algorithms: https://figshare.com/articles/Notes_on_pairwise_alignment_with_dynamic_programming/5223973
* Parasail: SIMD-accelerated pairwise alignment algorithms https://github.com/jeffdaily/parasail
* Greedy algorithm: http://pipmaker.bx.psu.edu/dist/greedy.pdf
* There is the possibility to apply techniques from short-read aligners to pairwise alignment, e.g. seed-and-chain. Some useful algorithms are available in Modular Aligner: https://github.com/ITBE-Lab/MA
* [WFA](https://github.com/smarco/WFA): `O(ns)` algorithm (i.e. very fast for similar sequences)
    * WFA-inspired implementation of LV https://github.com/lh3/lv89
* Block aligner (WIP): https://github.com/Daniel-Liu-c0deb0t/block-aligner
* https://github.com/CMU-SAFARI/Scrooge
* GPU-accelerated: https://github.com/fkallen/Accelign

### Local

* Nucl2vec: https://github.com/prakharg24/Nucl2vec
    * local alignment only - could it be adapted to semi-global?
    * Can re-implement in SpaCy? https://spacy.io/usage/vectors-similarity
* Accelerated version of Smith-Waterman algorithm implemented in the Genomics Kernel Library: https://github.com/Intel-HLS/GKL
* GPU-accelerated Smith-Waterman for nucleotides or proteins: https://github.com/mgawan/GPU-BSW

### Semi-global

[Semi-global alignment](
http://www.bioinf.uni-freiburg.de/Lehre/Courses/2013_SS/V_Bioinformatik_1/lecture4.pdf) is global alignment with an optional gap at the end(s).

* SeqPurge algorithm: insert match algorithm that performs thresholded exhaustive
   comparison to minimize probability of incorrect alignment. Relies on the fact that
   overlapping reads share alleles and indels (i.e. no gaps are required) (in C++).
   https://github.com/imgag/ngs-bits/tree/master/src/SeqPurge.
   * Speed up sequence comparison:
      * Between adapters and overhangs: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080745/pdf/btu177.pdf
      * Between reads: http://bioinformatics.oxfordjournals.org/content/30/14/2000.abstract
* Skewer algorithm: bit-masked k-difference matching (in C++): https://github.com/relipmoc/skewer
* Quality-aware overlap alignment (in Haskell): https://hackage.haskell.org/package/bio-0.5.3/docs/Bio-Alignment-QAlign.html
* FOGSAA, modified for semi-global alignment.
   * http://www.nature.com/articles/srep01746
   * http://www.isical.ac.in/~bioinfo_miu/FOGSAA.7z
* Adaptive banded alignment:
   * https://github.com/ocxtal/libgaba
   * https://github.com/ocxtal/adaptivebandbench
* Nepal: https://github.com/yamada-kd/nepal
* BGSA: provides accelerated implementation of the BitPAl algorithm https://github.com/sdu-hpcl/BGSA
* MEM-align: https://sites.google.com/site/memalignv1/
* Atria: byte-based matching algorithm https://github.com/cihga39871/Atria

### Global

* NW implementations with affine and linear gap penalties: https://github.com/AYahi/recNW
* A*PW https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner
* Parameter-free, context-specific scoring (optimized for long repeats): https://github.com/seryrzu/unialigner
* https://github.com/ruanjue/bsalign
* https://github.com/bxskdh/TSTA

### Barcode matching

* Pheniqs demultiplxer
   * Phred-adjusted ML for error probability: https://github.com/biosails/pheniqs/blob/master/docs/2.0/glossary.md

### Edit Distance

* https://github.com/maxdoblas/QuickEd

## Other

* Test whether different sequence encodings might enable faster alignment
https://github.com/hammerlab/kerseq/blob/master/kerseq/sequence_encoding.py
* Paper on relationship between alignment score and probability models: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz576/5536873
