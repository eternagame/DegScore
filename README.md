
Usage:

Clone `bpRNA` from https://github.com/hendrixlab/bpRNA.

add to arnie file:

`bprna: /path/to/bpRNA`


Example syntax:

```
from DegScore import DegScore
seq='GGGGAAAACCCC'
struct = '((((....))))'
mdl = DegScore(seq, struct)
print(mdl.counts)
print(mdl.score)
print(mdl.weights)
mdl.set_weights({'H': 1, 'E': 1, 'S': 1, 'I': 1, 'B': 1, 'M': 1})
print(mdl.score)
```

Output:
```
{'H': 4, 'E': 0, 'S': 8, 'I': 0, 'B': 0, 'M': 0} # mdl.counts
2.8 # mdl.score
{'H': 0.7, 'E': 1.0, 'S': 0.0, 'I': 0.2, 'B': 0.8, 'M': 0.6} # current mdl.weights
12 #updated mdl.score
```
