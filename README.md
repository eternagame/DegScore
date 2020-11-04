
Usage:

Clone `bpRNA` from https://github.com/hendrixlab/bpRNA.

add to arnie file:

`bprna: /path/to/bpRNA`

\[Update for Nov 2020 OpenVaccine work\]: DegScore class uses Arnie to calculate an MFE structure that is then used with bpRNA to compute loop type string.
DEFAULT setting is LinearFold-E (LinearFold with EternaFold) to mimic DegScore 2.1 calculation.

Example syntax:

```
from DegScore import DegScore

sequence='GGGGAAACCCC'

mdl = DegScore(sequence)
print('bpRNA string:')
print(mdl.bprna_string)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)
```

Output:
```
bpRNA string:
SSSSHHHSSSS
Predicted degscore per nucleotide:
[0.173 0.129 0.277 0.388 0.414 0.434 0.445 0.233 0.104 0.122 0.702]
Total degscore:
3.421000000000001
```

To change from these default LinearFold-E settings for MFE calculation to other Arnie settings (say, ViennaRNA):
```
mdl = DegScore(sequence, package='vienna', linear=False)
```
