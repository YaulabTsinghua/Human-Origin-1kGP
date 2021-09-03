Welcome to our project to explore the origin of human!

We propose a novel method named Divided Natural Vector baesd on the published Natural Vector method. The latter has been successfully applied on many datasets.

We strongly recommended to use `k=4` for human genomes analysis. The determination of k=4 is based on a dataset including 55 samples, as a subset of the dataset used in [1]. We tested the methods on the same dataset consisting of 55 samples. These samples include:
- 46 modern individuals
- 1 reference human genome
- 6 Neanderthals
-  2 other primates

The validation criterion was that the good approach had to distinguish between the following groups:
- human and other primates
- Neanderthals and modern individuals
- 15 modern African and 31 non-African individuals

The UPGMA trees for the natural vector method, the divided natural vector (k=4) method and MUSCLE are shown in Extended Data Figure 3(a), 3(b) and 3(c), respectively (can be found in `./tree`). All three approaches can distinguish between human and other primates, while the natural vector method fails to. Both the divided natural vector (k=4) method and MUSCLE get satisfactory results, but MUSCLE takes about 1 hour 18 minutes, while DNV takes only seconds. Therefore, we have proved that the divided natural vector method extracts the information hidden in the original sequences more accurately than the traditional Natural Vector method and is much more time-efficient compared to the alignment approach.

Due to the limitation of storage on Github, we add a toy data in `./data/` folder. You can find the recovered sequences for 2504 individiduals named as `xxIDxx.fa`, and with two vectors: the traditional natural vector and divided natural vector (k=4) in separate folders, both named as `out_xxIDxx.fa`. The computation is done in Perl with higher time efficiency compared to MATLAB.

The raw data of phylogenetic tree is stored in `./tree_files`, such as the nwk files and the figures. We will also upload all the trees (both in the main content and supplementary information) in this folder too.

Please contact yau@uic.edu if you have any questions about this project!


[1] Briggs A.W., Good J.M., Green R.E., et al. Targeted retrieval and analysis of five Neandertal mtDNA genomes. 2009. Science, 325:318-321.
