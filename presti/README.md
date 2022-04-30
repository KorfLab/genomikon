README for the Prestidigitizer
==============================

Number stream classifier.

## Background ##

The Prestidigitizer, `presti`, is a classifier for streams of numbers. The
original intent of the program is for peak finding in ChIP-seq-like data.
Assume you are given a stream of numbers like the following:

```
1.0
1.1
1.1
1.0
9.9
9.8
9.9
9.7
1.0
1.0
1.2
1.0
```

There are two different classes of coverage, one whose value is around 1 and
another whose value is around 10. Given that there appear to be two classes,
you would tell `presti` to decode the stream as two classes: "1" and "10".
After running `presti`, it will tell you there are three regions:

```
# beg   end     class
1       4       1
5       8       10
9       12      1
```

## Tuning ##

* `-x` controls splitting vs. merging
* `-y` controls the scale of numbers
* `-z` limits the maximum value

The `-x` parameter is the cost to switch from one state to another. The higher
the value of `-x` the harder it is for the class to switch. If you are finding
multiple peaks that are very close together and want to merge them into single
peaks, increase `-x`.

The `-y` parameter is used to scale up or down numbers. This simply multiplies
all input values by `-y`. This might be useful if the input numbers are very
high because `presti` is limited to a maximum of 999.

The `-z` parameter limits the maximum size of a number, which is usually 999.
Any value above `-z` gets converted to `-z`.


## Model ##

The `presty` hidden Markov model (HMM) is fully connected among all states. The
transition probabilities for staying in the same state all have the same
relatively high value and the transition probabilities for switching states all
have the same relatively low value. The low transition probabilities correspond
to whatever the `-x` value is set to. For example, `-x 5` corresponds to about
0.0067 (1/e^5). Theoretically, in a 2-state HMM, the self-returning transition
should therefore be 1 - 1/e^5, but `presti` uses 1.0 for all self-returning
transitions regardless of the number of states or the value of `-x`. If that
sounds offensive, you can call the `presti` model a conditional random field
(CRF) instead of an HMM.

The emission probabilities within each state are drawn from a Poisson
distribution. Consequently, a class with an expected value of 1 will generate
1s much more often than 10s, and a class with an expected value of 10 will
generate 10s much more often than 1s. Since the emissions are drawn from a
Poisson distribution, they naturally sum to 1.0, so at least this isn't
offensive (unless you consider that all of the calculations are performed in
log-space).

The `presti` uses the Viterbi algorithm for decoding the input file. It
therefore produces the maximum likelihood path, but doesn't give you any
information about the probability of that path.
