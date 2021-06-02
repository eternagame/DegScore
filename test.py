from DegScore import DegScore

sequence='GGGUUUCCC'
structure='(((...)))'
mdl = DegScore(sequence, structure=structure)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)

hkws_result_02jun2021 = [0.17,  0.217, 0.511, 0.641, 0.714, 0.49,  0.347, 0.15,  0.643]

for i in range(len(sequence)):
	assert abs(mdl.degscore_by_position[i] - hkws_result_02jun2021[i]) < 1e-5
