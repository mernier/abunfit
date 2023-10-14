#!/usr/bin/python

import os
import sys
import math
import numpy


def main():
	file=sys.argv[1]
	outputfile=sys.argv[2]


	f=open(file,'r')
	lines=f.readlines()
	number_models = len(lines[0].split(" ")) - 2


	# # Read data...

	chi2=numpy.arange(len(lines))
	chi2=chi2*0.0
	alpha=numpy.arange(len(lines))
	alpha=alpha*0.0
	model1 = [None]*len(lines)
	if number_models > 1:
		model2 = [None]*len(lines)
	if number_models > 2:
		model3 = [None]*len(lines)
	if number_models > 3:	
		model4 = [None]*len(lines)
	if number_models > 4:
		sys.exit("Error: too many colums in the input file (4 models max.)")


	i=0
	for line in lines:
		columns=line.split(" ")
		chi2[i]=float(columns[0])
		alpha[i]=float(columns[1])
		model1[i]=columns[2]
		if number_models > 1:
			model2[i]=columns[3]
		if number_models > 2:
			model3[i]=columns[4]
		if number_models > 3:
			model4[i]=columns[5]
		i=i+1






	# # Initializing variables...

	chi2_sorted=chi2*0.0
	alpha_sorted=alpha*0.0
	model1_sorted=[None]*len(lines)
	if number_models > 1:
		model2_sorted=[None]*len(lines)
	if number_models > 2:
		model3_sorted=[None]*len(lines)
	if number_models > 3:
		model4_sorted=[None]*len(lines)





	# # Sorting process

	# First loop
	oldchi=1000000000000000
	oldalpha=0.0
	oldmodel1='None'
	if number_models > 1:
		oldmodel2='None'
	if number_models > 2:
		oldmodel3='None'
	if number_models > 3:
		oldmodel4='None'
	for i in range(len(chi2)):
			if chi2[i] < oldchi:
				oldchi=chi2[i]
				oldalpha=alpha[i]
				oldmodel1=model1[i]
				if number_models > 1:
					oldmodel2=model2[i]
				if number_models > 2:
					oldmodel3=model3[i]
				if number_models > 3:
					oldmodel4=model4[i]


	chi2_sorted[0]=oldchi
	alpha_sorted[0]=oldalpha
	model1_sorted[0]=oldmodel1
	if number_models > 1:
		model2_sorted[0]=oldmodel2
	if number_models > 2:
		model3_sorted[0]=oldmodel3
	if number_models > 3:
		model4_sorted[0]=oldmodel4





	# Other loops
	for j in range(1,len(chi2)):
		oldchi=1000000000000000000
		oldalpha=0.0
		oldmodel1='None'
		if number_models > 1:
			oldmodel2='None'
		if number_models > 2:
			oldmodel3='None'
		if number_models > 3:
			oldmodel4='None'

		for i in range(len(chi2)):
				if chi2_sorted[j-1] < chi2[i] < oldchi:
					oldchi=chi2[i]
					oldalpha=alpha[i]
					oldmodel1=model1[i]
					if number_models > 1:
						oldmodel2=model2[i]
					if number_models > 2:
						oldmodel3=model3[i]
					if number_models > 3:
						oldmodel4=model4[i]
		chi2_sorted[j]=oldchi
		alpha_sorted[j]=oldalpha
		model1_sorted[j]=oldmodel1
		if number_models > 1:
			model2_sorted[j]=oldmodel2
		if number_models > 2:
			model3_sorted[j]=oldmodel3
		if number_models > 3:
			model4_sorted[j]=oldmodel4		



	chi2_sorted = chi2_sorted/100000.0	#...because in input, chi2 was already multiplied by 1000...






	# # Printing...

	fileout = open(outputfile, 'w+')

	for i in range(len(lines)):

		if number_models == 1:
			fileout.write("%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i]))
		if number_models == 2:
			fileout.write("%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i]))
		if number_models == 3:
			fileout.write("%s	%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i], model3_sorted[i]))
		if number_models == 4:
			fileout.write("%s	%s	%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i], model3_sorted[i], model4_sorted[i]))
  


	fileout.close()
	f.close()





main()
