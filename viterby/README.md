viterby
=======

A demonstration program of how to write a viterbi decoder for an HMM specified
in JSON format.

HMM file format
---------------

The file format is JSON, but it must be crafted very specifically because I'm 
not using an actual JSON parser. The overall structure of an HMM is simple as 
it is just a collection of states.

```
{
	"name": "exon-intron",
	"states": 2,
	"state": [
		{state_object}
		{state_object}
	]
}
```

State objects are a bit more complicated. The number of outgoing transitions is 
under `transitions`. The total number of emission probabilities is under 
`emissions`, which should be 4**n, where _n_ is the order of the Markov model. 
For explicit durations, the length of the defined region and the geometric tail 
are under `durations` and `tail`. The values for transitions, emissions, and 
durations follow, as shown below.


```
{
	"name": "exon",
	"init": 0.5,
	"term": 0.5,
	"transitions": 2,
	"emissions": 4,
	"durations": 0,
	"tail": 0,
	"trans": {
		"exon": 0.99,
		"intron": 0.01
	}
	"emits": [0.2, 0.3, 0.3, 0.2],
	"durations": [],
}
```
