#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def Count(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum


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

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity', 'polarity','charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in range(1, len(groups) + 1):
			header.append(p + '.G' + str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for p in property:
			c1 = Count(group1[p], sequence) / len(sequence)
			c2 = Count(group2[p], sequence) / len(sequence)
			c3 = 1 - c1 - c2
			code = code + [c1, c2, c3]
		encodings.append(code)
	return encodings
