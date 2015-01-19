HHCompare -- code for HMM based paralogues analysis
=========================================================

HHCompare is a pipeline for HMM-HMM comparison based hierarchial clustering 
and analysis of potential paralogues in sequence set.
It is based on premise that paralogues sequences have very high similarity
(below certain e-value cutoff for HMM vs HMM alignment). 

code generates HMM models from input (set of sequences), compares them
all vs all, groups similar ones and then repeats; simplified workflow is as following:

1. generate HMMs from sequences or groups
2. do HMM-HMM comparison between HMMs generated under 1)
3. evalute pairs for similarity, merge pair(s) with similarity below CUTOFF in group(s)
4. repeat 1) if any merge was done, finish run if no merge was done
5. parse results and generate hierarchial trees of groups 

Notes: 
- input: multiple protein sequences in FASTA format; not tested for nucleic acid sequences
- output: set of trees respresenting potential paralogue clusters and unclustered sequences
- note that code DOES NOT provide graphical representation of results, it generates newick outputs 
for detected groups of clustered sequences; various tools can convert it into tree-like graphics (example: http://etetoolkit.org/treeview/)

Dependencies: 
- RBioTools.py : various bioinformatics functions; should be in same folder as code
- HH-suite2.0, installed locally; make sure paths are set property when running
    download at ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/, follow appropriate instructions for install
- ClustalW2, installed locally;make sure paths are set property when running
    - download at ftp://ftp.ebi.ac.uk/pub/software/clustalw2/
- python2.7

More notes: 
- last update: 19/01/2015
- if python is not in /usr/bin/python2.7, code will not self-execute; run it as <python> <codename> or edit 1st line
- in order to avoid entering paths each time code is run, change lines 329-331 default= to appropriate default
- tested under Ubuntu 12.04.5 LTS
- to annotate results, use HHCompareAnnotator.py (note: script is in early prototype version, and might need tweaking)

Example run: 
- ./HHCompare.py -I ./test.fa -O __test_out --HHMAKE <path> --HHALIGN <path> --CLUSTALW2 <path> 