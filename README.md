genomikon
=========

Genomic sequence analysis library primarily by Ian Korf.

## Philosophy ##

+ Simple & Native
+ Noob-friendly
+ OOP-ish

### Simple & Native

There aren't a whole bunch of typedefs that rename built-in types. DNA sequences
are just strings. There is no Boolean type.

### Noob-friendly

The code is generally not overly sophisticated or brief. Verbosity in code is
encouraged. Verbosity in comments is not. Let the compiler do the optimizations
and let the coder read the code.

## OOP-ish ##

Every complex data structure, aka _object_, is a pointer to a struct. There is
always a constructor function (generally named with `new` or `read`) that
creates and initializes the object and a `free` function that destroys the
object.

### No side-effects

There are no function side-effects. Everything going into a function is generally `const`. Everything coming back from a function is either a built-in type or an object (pointer to a struct).

## Installation ##

Just `make` from inside the `genomikon` directory to compile the library.

	make

If you want to build the demo apps, you can `make demo` or change to the `demo`
directory and `make` from there.

The `testing` program is used to check for memory leaks and such. It's not
necessary, but you can `make test` if you want.

## Core Functions ##

### Vectors

Vectors - dynamic arrays

+ vec - void * vector
+ ivec - integer vector
+ fvec - double precisions floating point vector
+ tvec - text vector

### Maps

Maps - dictionaries, hashes

+ map - void * map
+ tmap - text map

### FASTA I/O

### GFF I/O

### Sequence Models

+ pwm - position weight matrix
+ mm - Markov model
+ len - length model

### Hidden Markov Models

HMMs are defined by their states.

### Command line processing



## Demo Programs ##

There are a few demonstration programs in the `demo` directory.

+ dusty - low complexity filter (CLI, fasta files)
+ geney - gene scoring program (CLI, fasta, gff, models)
+ smithy - Smith-Waterman alignment
+ viterby - gene prediction

## To Do ##

- add HMMs
