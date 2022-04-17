#!/usr/bin/env python
#_*_coding:utf-8_*_

import re , math

def Count(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum
def Count1(aaSet, sequence):
	number = 0
	for aa in sequence:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(sequence)):
			if sequence[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(sequence) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code

def CTDC(fastas, **kw):
	group1 = {
		'hydrophobicity': 'RKEDQN',
		'polarity':        'LIFWCMVY',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity': 'GASTPHY',
		'polarity':        'PATGS',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity': 'CLVIMFW',
		'polarity':        'HQRKNED',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}
	AAGroup = {
		'g1': 'DE',
		'g2': 'HRK',
		'g3': 'CGNQSTY',
		'g4': 'AFILMPVW'
		
	}

	myGroups = sorted(AAGroup.keys())

	AADict = {}
	for g in myGroups:
		for aa in AAGroup[g]:
			AADict[aa] = g

	features = [f1 + '.'+ f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity', 'polarity','charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in range(1, len(groups) + 1):
			header.append(p + '.G' + str(g))
		for tr in ('Tr1', 'Tr2', 'Tr3'):
			header.append(p + '.' + tr)
		for g in ('1', '2', '3'):
			for d in ['0', '25', '50', '75', '100']:
				header.append(p + '.G' + g + '.res@' + d)
	encodings.append(header)
	
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
		for p in property:
			c1 = Count(group1[p], sequence) / len(sequence)
			c2 = Count(group2[p], sequence) / len(sequence)
			c3 = 1 - c1 - c2
			c1221, c1331, c2332 = 0, 0, 0
			for pair in aaPair:
				if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
					c1221 = c1221 + 1
					continue
				if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
					c1331 = c1331 + 1
					continue
				if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
					c2332 = c2332 + 1
			code = code + [c1, c2, c3]+[c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)] + Count1(group1[p], sequence) + Count1(group2[p], sequence) + Count1(group3[p], sequence)
		encodings.append(code)
	return encodings

	
