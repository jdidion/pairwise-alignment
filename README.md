# Pairwise alignment algorithms

## Semi-global

[Semi-global alignment](
http://www.bioinf.uni-freiburg.de/Lehre/Courses/2013_SS/V_Bioinformatik_1/lecture4.pdf) is global alignment with an optional gap at the end(s).

1. SeqPurge algorithm: insert match algorithm that performs thresholded exhaustive
   comparison to minimize probability of incorrect alignment. Relies on the fact that
   overlapping reads share alleles and indels (i.e. no gaps are required) (in C++).
   https://github.com/imgag/ngs-bits/tree/master/src/SeqPurge.
   * Speed up sequence comparison:
     * Between adapters and overhangs:
       * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080745/pdf/btu177.pdf
     * Between reads:
       * http://bioinformatics.oxfordjournals.org/content/30/14/2000.abstract
2. Skewer algorithm: bit-masked k-difference matching (in C++): https://github.com/relipmoc/skewer
3. Quality-aware overlap alignment (in Haskell): https://hackage.haskell.org/package/bio-0.5.3/docs/Bio-Alignment-QAlign.html
4. FOGSAA, modified for semi-global alignment.
   * http://www.nature.com/articles/srep01746
   * http://www.isical.ac.in/~bioinfo_miu/FOGSAA.7z
5. EDLIB: edit distance-based alignment: https://github.com/Martinsos/edlib
6. Phred-adjusted ML for error probability:
https://biosails.github.io/pheniqs/glossary.html#phred_adjusted_maximum_likelihood_decoding
7. Adaptive banded alignment: https://github.com/ocxtal/libgaba
8. Nepal: https://github.com/yamada-kd/nepal
9. The SeqAn C++ library implements several alignment algorithms:
http://www.sciencedirect.com/science/article/pii/S0168165617315420
10. Alignment-free tools for pairwise sequence comparison:
* http://www.combio.pl/alfree/tools/
* http://www.combio.pl/alfree
* http://bioinformatics.org.au/tools/decaf+py/
11. https://github.com/sdu-hpcl/BGSA
12. Train a NN to approximate the pairwise alignment distance
https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty887/5140215?redirectedFrom=fulltext
13. Nucl2vec: https://github.com/prakharg24/Nucl2vec (local alignment only - could
it be adapted to semi-global?)
* Can re-implement in SpaCy? https://spacy.io/usage/vectors-similarity
14. Accelerated version of Smith-Waterman algorithm implemented in the Genomics
Kernel Library: https://github.com/Intel-HLS/GKL
15. MEM-align: https://sites.google.com/site/memalignv1/

## Other

* Test whether different sequence encodings might enable faster alignment
https://github.com/hammerlab/kerseq/blob/master/kerseq/sequence_encoding.py
