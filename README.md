
Usage:

Clone `bpRNA` from https://github.com/hendrixlab/bpRNA.

add to arnie file:

`bprna: /path/to/bpRNA`


Example syntax:

```
from DegScore import DegScore

sequence='GGGGAAACCCC'

mdl = DegScore(sequence)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)
```

Output:
```
Predicted degscore per nucleotide:
[0.173 0.129 0.277 0.388 0.414 0.434 0.445 0.233 0.104 0.122 0.702]
Total degscore:
3.421000000000001
```
