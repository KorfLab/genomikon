genomikon
=========

Genomic sequence analysis library and applications. Not intended to be used for
large projects or lab projects. More like Ian's private/public C sandbox.


Installation
------------

Building the software is as simple as:

	make

This will build the genomikon library and the applications below. To perform
some functional tests, you can:

	make test

The `testing` program in the main directory is used to check for memory leaks
and is not part of the functional tests.


Applications
------------

There are several applications built with the genomikon library. Each of these
has its own directory.

+ dusty - a complexity filter for nucleotide sequences
+ isoformer - programs for generating alternative splicing isoforms
+ hmmstar - a viterbi decoder for nucleotide sequences and k-mer models
+ presti - an HMM for classifying integers (e.g. read coverage)
+ smithy - a simple implementation of Smith-Waterman
+ viterby - a simple Viterbi decoder for HMMs specified in JSON
+ wordy - a word game that mixes boggle and genetic algorithms


Author
------

Ian Korf
