genomikon
=========

Genomic sequence analysis

## Installation ##

Just `make` from inside the `genomikon` directory to compile the library.

	make

If you want to build the demo apps, you can `make demo` or change to the `demo`
directory and `make` from there.

The `testing` program is used to check for memory leaks and such. It's not
necessary, but you can `make test` if you want.

## Core Functions ##

## Demo Programs ##

There are a few demonstration programs in the `demo` directory.

+ dusty - low complexity filter
+ geney - gene scoring program (may change to gene prediction)
+ smithy - Smith-Waterman alignment

## To Do ##

- Change all file i/o to pipe and readline
- Move denada into genomikon library
