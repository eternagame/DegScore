from DegScore import DegScore

sequence='GGGGAAACCCC'

mdl = DegScore(sequence)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)

# Test output matches hkws notebook calculation of DegScore2.1

hkws_result_from_notebook = [0.173, 0.129, 0.277, 0.388, 0.414,
       0.434, 0.445, 0.233, 0.104, 0.122, 0.702]

for i in range(len(sequence)):
	assert abs(mdl.degscore_by_position[i] - hkws_result_from_notebook[i]) < 1e-5
