viterby
=======

The HMM file format is JSON, but it is crafted in a line-specific manner in 
case the language doesn't have a built-in JSON parser. The overall structure of 
the HMM is some meta-data and then an array of states (see below).

```
{
	"name": "exon-intron",
	"author": "Ian",
	"version": 1.0,
	"comments": "this is for testing",
	"states": 2,
	"state": [
		{state_object}
		{state_object}
	]
}
```

State objects are a bit more complicated. Forward pointing transitions are 
stored in the `transitions` object. The value for `emissions` should be 4**n 
where n is the order of the Markov model, and the specific values are stored in 
a related array.


```
{
	"name": "exon",
	"init": 0.5,
	"term": 0.5,
	"transitions": 2,
		"transition": {
		"exon": 0.99,
		"intron": 0.01
	}
	"emissions": 4,
	"emission": [0.2, 0.3, 0.3, 0.2],
}
```

States with explicit durations (not common) have tags for `durations`, 
`duration`, and `tail`. These can be given uninteresting values (as shown) or 
omitted entirely.

```
{
	"name": "exon",
	"init": 0.5,
	"term": 0.5,
	"transitions": 2,
		"transition": {
		"exon": 0.99,
		"intron": 0.01
	}
	"emissions": 4,
	"emission": [0.2, 0.3, 0.3, 0.2],
	"durations": 0,
	"duration": false,
	"tail": false,
}

