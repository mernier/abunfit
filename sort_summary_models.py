#!/usr/bin/python

import os
import sys
import math
import numpy


def main():
	file=sys.argv[1]


	f=open(file,'r')
	lines=f.readlines()



	# # Read data...

	chi2=numpy.arange(len(lines))
	chi2=chi2*0
	alpha=numpy.arange(len(lines))
	alpha=alpha*0.0
	model1 = [None]*len(lines)
	model2 = [None]*len(lines)
	model3 = [None]*len(lines)
	model4 = [None]*len(lines)

	i=0
	for line in lines:
		columns=line.split(" ")
		chi2[i]=float(columns[0])
		alpha[i]=float(columns[1])
		model1[i]=columns[2]
		if len(columns)>3:
			model2[i]=columns[3]
		if len(columns)>4:
			model3[i]=columns[4]
		if len(columns)>5:
			model4[i]=columns[5]
		i=i+1






	# # Initializing variables...

	chi2_sorted=chi2*0
	alpha_sorted=alpha*0.0
	model1_sorted=[None]*len(lines)
	model2_sorted=[None]*len(lines)
	model3_sorted=[None]*len(lines)
	model4_sorted=[None]*len(lines)




	# # Sorting process

	# First loop
	oldchi=1000000000000000
	oldalpha=0.0
	oldmodel1='None'
	oldmodel2='None'
	oldmodel3='None'
	oldmodel4='None'
	for i in range(len(chi2)):
			if chi2[i] < oldchi:
				oldchi=chi2[i]
				oldalpha=alpha[i]
				oldmodel1=model1[i]
				if len(columns)>3:
					oldmodel2=model2[i]
				if len(columns)>4:
					oldmodel3=model3[i]
				if len(columns)>5:
					oldmodel4=model4[i]

	chi2_sorted[0]=oldchi
	alpha_sorted[0]=oldalpha
	model1_sorted[0]=oldmodel1
	if len(columns)>3:
		model2_sorted[0]=oldmodel2

	if len(columns)>4:
		model3_sorted[0]=oldmodel3

	if len(columns)>5:
		model4_sorted[0]=oldmodel4





	# Other loops
	for j in range(1,len(chi2)):
		oldchi=1000000000000000000
		oldalpha=0.0
		oldmodel1='None'
		oldmodel2='None'
		oldmodel3='None'
		oldmodel4='None'
		for i in range(len(chi2)):
				if chi2_sorted[j-1] < chi2[i] < oldchi:
					oldchi=chi2[i]
					oldalpha=alpha[i]
					oldmodel1=model1[i]
					if len(columns)>3:
						oldmodel2=model2[i]
					if len(columns)>4:
						oldmodel3=model3[i]
					if len(columns)>5:
						oldmodel4=model4[i]
		chi2_sorted[j]=oldchi
		alpha_sorted[j]=oldalpha
		model1_sorted[j]=oldmodel1
		if len(columns)>3:
			model2_sorted[j]=oldmodel2
		if len(columns)>4:
			model3_sorted[j]=oldmodel3
		if len(columns)>5:
			model4_sorted[j]=oldmodel4


	chi2_sorted = chi2_sorted/100000.0	#...because in input, chi2 was already multiplied by 1000...






	# # Printing...

	fileout = open("summary_models_sorted.txt", 'w+')

	if len(columns)==3:
		for i in range(len(lines)):
			fileout.write("%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i]))

	if len(columns)==4:
		for i in range(len(lines)):
			fileout.write("%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i]))
	
	if len(columns)==5:
		for i in range(len(lines)):
			fileout.write("%s	%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i], model3_sorted[i]))

	if len(columns)==6:
		for i in range(len(lines)):
			fileout.write("%s	%s	%s	%s	%s	%s" % (chi2_sorted[i], alpha_sorted[i], model1_sorted[i], model2_sorted[i], model3_sorted[i], model4_sorted[i]))

  


	fileout.close()
	f.close()





main()
