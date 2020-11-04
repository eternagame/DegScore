from DegScore import DegScore

sequence='GGGGAAACCCC'

mdl = DegScore(sequence)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)

# Test output matches hkws notebook calculation of DegScore2.1

hkws_result_from_notebook = [0.17303802, 0.12903802, 0.27703802, 0.38803802, 0.41403802,
       0.43403802, 0.44503802, 0.23303802, 0.10403802, 0.12203802, 0.70203802]

for i in range(len(sequence)):
	assert abs(mdl.degscore_by_position[i] - hkws_result_from_notebook[i]) < 1e-5
