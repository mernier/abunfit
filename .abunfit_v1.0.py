#!/usr/bin/python

import os
import sys
import math
import numpy
import scipy.optimize as optimization


def main():

	if len(sys.argv) < 4:
		sys.exit("Error: Uncorrect arguments.\n\n \
		PROGRAM: abunfit.py\n\n \
		SYNOPSIS: ./abunfit.py [input_file] [number_of_models] [alpha] [SNe_model_1] [SNe_model_2] ... \n\n \
		DESCRIPTION: Fits a set of intra-cluster abundances with a combination of \n \
		supernovae (SNe) yields models. Several models (either SNcc or SNIa) are available \n \
		and can be updated or created. \n\n \
		ARGUMENTS: \n\n \
		        [input_file] -- Text input file of the set of measured ICM abundances \n \
		                        and that contains 3 columns: (1) the element number (Z), \n \
		                        (2) the measured abundance, and (3) the associated \n \
		                        (symmetrical) uncertainties.\n\n \
		        [number_of_models] -- Integer (from 1 to 4) indicating the number of \n \ 			                              models to be fitted simultaneously (min: 2, max: 4).\n \
		        [alpha] -- Slope index of the assumed initial mass function (IMF). \n \
		                   Ex: for Salpeter IMF, alpha = -2.35.\n \
		        [model<1/2/...>] -- Name of the model. Currently available models: \n \
		                            ==== SNIa ====\n \
		                            == Iwamoto et al. (1999)\n \
		                            - Iwamoto_W7\n \
		                            - Iwamoto_W70\n \
		                            - Iwamoto_WDD1\n \
		                            - Iwamoto_WDD2\n \
		                            - Iwamoto_WDD3\n \
		                            - Iwamoto_CDD1\n \
		                            - Iwamoto_CDD2\n \
		                            == Badenes et al. (2006)\n \
		                               (Based on Tycho's best models...)\n \
		                            - Badenes_DDTa\n \
		                            - Badenes_DDTb\n \
		                            - Badenes_DDTc\n \
		                            - Badenes_DDTd\n \
		                            - Badenes_DDTe\n \
		                            - Badenes_DDTf\n \
		                            == Maeda et al. (2010)\n \
		                               (2-D explosion models...)\n \
		                            - Maeda_C-DEF\n \
		                            - Maeda_C-DDT\n \
		                            - Maeda_O-DDT\n \
		                            ==== Ca-rich gap transient SNe ====\n \
		                            == Waldman et al. (2011)\n \
		                            - Waldman_CO.45HE.2\n \
		                            - Waldman_CO.55HE.2\n \
		                            - Waldman_CO.5HE.15\n \
		                            - Waldman_CO.5HE.2\n \
		                            - Waldman_CO.5HE.2C.3\n \
		                            - Waldman_CO.5HE.2N.02\n \
		                            - Waldman_CO.5HE.3\n \
		                            - Waldman_CO.6HE.2\n \
		                            ==== SNcc ====\n \
		                            == Nomoto et al. (2006)\n \
		                               (One initial metallicity per model...)\n \
		                            - 0\n \
		                            - 0.001\n \
		                            - 0.004\n \
		                            - 0.02\n\n \
		EXAMPLE: ./abunfit.py average_abun_clusters.out 2 -2.35 Iwamoto_WDD2 0.02\n\n \
		VERSION: 1.1 (January 2016)\n \
		         1.0 (July 2015)\n\n \
		AUTHOR: Francois Mernier\n\n")




	inputfile=sys.argv[1]
	number_of_models=sys.argv[2]
	number_of_models = int(number_of_models)
	alpha=sys.argv[3]
	input_mod_1=sys.argv[4]
	if number_of_models>1:
		input_mod_2=sys.argv[5]
	if number_of_models>2:
		input_mod_3=sys.argv[6]
	if number_of_models>3:
		input_mod_4=sys.argv[7]
	if number_of_models>4 or number_of_models<1:
		sys.exit("Error: This program can only fit 1,2,3 or 4 models simultaneously... Sorry :(")




	# Read cluster measurements data.
	f1 = file(inputfile)
	cl_abun = numpy.loadtxt(f1)

	elements = cl_abun[:,0]
	data = cl_abun[:,1]
	errors = cl_abun[:,2]

	# IMF assumption: Salpeter
	#alpha = -2.35
	aplha = float(alpha)



	# Give to each considered element its respective name.
	elem_name = [None]*len(elements)
	for i in range(len(elements)):
		if int(elements[i]) == 1:
			elem_name[i] = 'H'
		if int(elements[i]) == 2:
			elem_name[i] = 'He'
		if int(elements[i]) == 3:
			elem_name[i] = 'Li'
		if int(elements[i]) == 4:
			elem_name[i] = 'Be'
		if int(elements[i]) == 5:
			elem_name[i] = 'B'
		if int(elements[i]) == 6:
			elem_name[i] = 'C'
		if int(elements[i]) == 7:
			elem_name[i] = 'N'
		if int(elements[i]) == 8:
			elem_name[i] = 'O'
		if int(elements[i]) == 9:
			elem_name[i] = 'F'
		if int(elements[i]) == 10:
			elem_name[i] = 'Ne'
		if int(elements[i]) == 11:
			elem_name[i] = 'Na'
		if int(elements[i]) == 12:
			elem_name[i] = 'Mg'
		if int(elements[i]) == 13:
			elem_name[i] = 'Al'
		if int(elements[i]) == 14:
			elem_name[i] = 'Si'
		if int(elements[i]) == 15:
			elem_name[i] = 'P'
		if int(elements[i]) == 16:
			elem_name[i] = 'S'
		if int(elements[i]) == 17:
			elem_name[i] = 'Cl'
		if int(elements[i]) == 18:
			elem_name[i] = 'Ar'
		if int(elements[i]) == 19:
			elem_name[i] = 'K'
		if int(elements[i]) == 20:
			elem_name[i] = 'Ca'
		if int(elements[i]) == 21:
			elem_name[i] = 'Sc'
		if int(elements[i]) == 22:
			elem_name[i] = 'Ti'
		if int(elements[i]) == 23:
			elem_name[i] = 'V'
		if int(elements[i]) == 24:
			elem_name[i] = 'Cr'
		if int(elements[i]) == 25:
			elem_name[i] = 'Mn'
		if int(elements[i]) == 26:
			elem_name[i] = 'Fe'
		if int(elements[i]) == 27:
			elem_name[i] = 'Co'
		if int(elements[i]) == 28:
			elem_name[i] = 'Ni'
		if int(elements[i]) == 29:
			elem_name[i] = 'Cu'
		if int(elements[i]) == 30:
			elem_name[i] = 'Zn'
		#else:
		#	print "Error: could not find the name of element", int(elements[i])






	# Read atomic masses. 
	# mass[i] = atomic mass of the ith element.
	f2 = file("atomic_masses.txt")
	atom_mass = numpy.loadtxt(f2)

	mass = numpy.arange(len(elements))
	mass = mass*0.0
	for i in range(len(elements)):
		mass[i] = atom_mass[elements[i]-1,1]





	# Read abundance tables (Lodders et al. 2009). 
	# lodders[i] = protosolar abun. of the ith element (rel. to Si==10^6).
	f3 = file("lodders09.txt")
	abun_table = numpy.loadtxt(f3)

	lodders = numpy.arange(len(elements))
	lodders = lodders*0.0
	for i in range(len(elements)):
		lodders[i] = abun_table[elements[i]-1,1]






	# Load SNIa and/or SNcc models...

	if (input_mod_1 == "Iwamoto_W7") or (input_mod_1 == "Iwamoto_W70") or \
		(input_mod_1 == "Iwamoto_WDD1") or (input_mod_1 == "Iwamoto_WDD2") or \
		(input_mod_1 == "Iwamoto_WDD3") or (input_mod_1 == "Iwamoto_CDD1") or \
		(input_mod_1 == "Iwamoto_CDD2") or \
		(input_mod_1 == "Badenes_DDTa") or (input_mod_1 == "Badenes_DDTb") or \
		(input_mod_1 == "Badenes_DDTc") or (input_mod_1 == "Badenes_DDTd") or \
		(input_mod_1 == "Badenes_DDTe") or (input_mod_1 == "Badenes_DDTf") or \
		(input_mod_1 == "Badenes_DDTa_undecayed") or (input_mod_1 == "Badenes_DDTb_undecayed") or \
		(input_mod_1 == "Badenes_DDTc_undecayed") or (input_mod_1 == "Badenes_DDTd_undecayed") or \
		(input_mod_1 == "Badenes_DDTe_undecayed") or (input_mod_1 == "Badenes_DDTf_undecayed") or \
		(input_mod_1 == "Maeda_C-DEF") or (input_mod_1 == "Maeda_C-DDT") or \
		(input_mod_1 == "Maeda_O-DDT") or \
		(input_mod_1 == "Waldman_CO.45HE.2") or (input_mod_1 == "Waldman_CO.55HE.2") or \
		(input_mod_1 == "Waldman_CO.5HE.15") or (input_mod_1 == "Waldman_CO.5HE.2") or \
		(input_mod_1 == "Waldman_CO.5HE.2C.3") or (input_mod_1 == "Waldman_CO.5HE.2N.02") or \
		(input_mod_1 == "Waldman_CO.5HE.3") or (input_mod_1 == "Waldman_CO.6HE.2"):
		model1 = sum_SNIa(input_mod_1,elements)
	elif (input_mod_1 == "0") or (input_mod_1 == "0.001") or \
		(input_mod_1 == "0.004") or (input_mod_1 == "0.02") or \
		(input_mod_1 == "0_undecayed") or (input_mod_1 == "0.001_undecayed") or \
		(input_mod_1 == "0.004_undecayed") or (input_mod_1 == "0.02_undecayed"):
		model1 = integrate_SNcc(input_mod_1,alpha,elem_name)
	else:
		print "Error: The model", input_mod_1, "does not exit, or could not be found."
		sys.exit()


	if number_of_models>1:
		if (input_mod_2 == "Iwamoto_W7") or (input_mod_2 == "Iwamoto_W70") or \
			(input_mod_2 == "Iwamoto_WDD1") or (input_mod_2 == "Iwamoto_WDD2") or \
			(input_mod_2 == "Iwamoto_WDD3") or (input_mod_2 == "Iwamoto_CDD1") or \
			(input_mod_2 == "Iwamoto_CDD2") or \
			(input_mod_2 == "Badenes_DDTa") or (input_mod_2 == "Badenes_DDTb") or \
			(input_mod_2 == "Badenes_DDTc") or (input_mod_2 == "Badenes_DDTd") or \
			(input_mod_2 == "Badenes_DDTe") or (input_mod_2 == "Badenes_DDTf") or \
			(input_mod_2 == "Badenes_DDTa_undecayed") or \
			(input_mod_2 == "Badenes_DDTb_undecayed") or \
			(input_mod_2 == "Badenes_DDTc_undecayed") or \
			(input_mod_2 == "Badenes_DDTd_undecayed") or \
			(input_mod_2 == "Badenes_DDTe_undecayed") or \
			(input_mod_2 == "Badenes_DDTf_undecayed") or \
			(input_mod_2 == "Maeda_C-DEF") or \
			(input_mod_2 == "Maeda_C-DDT") or \
			(input_mod_2 == "Maeda_O-DDT") or \
			(input_mod_2 == "Waldman_CO.45HE.2") or \
			(input_mod_2 == "Waldman_CO.55HE.2") or \
			(input_mod_2 == "Waldman_CO.5HE.15") or \
			(input_mod_2 == "Waldman_CO.5HE.2") or \
			(input_mod_2 == "Waldman_CO.5HE.2C.3") or \
			(input_mod_2 == "Waldman_CO.5HE.2N.02") or \
			(input_mod_2 == "Waldman_CO.5HE.3") or \
			(input_mod_2 == "Waldman_CO.6HE.2"):
			model2 = sum_SNIa(input_mod_2,elements)
		elif (input_mod_2 == "0") or (input_mod_2 == "0.001") or \
			(input_mod_2 == "0.004") or (input_mod_2 == "0.02") or \
			(input_mod_2 == "0_undecayed") or (input_mod_2 == "0.001_undecayed") or \
			(input_mod_2 == "0.004_undecayed") or (input_mod_2 == "0.02_undecayed"):
			model2 = integrate_SNcc(input_mod_2,alpha,elem_name)
		else:
			print "Error: The model", input_mod_2, "does not exit, or could not be found."
			sys.exit()


	if number_of_models>2:
		if (input_mod_3 == "Iwamoto_W7") or (input_mod_3 == "Iwamoto_W70") or \
			(input_mod_3 == "Iwamoto_WDD1") or (input_mod_3 == "Iwamoto_WDD2") or \
			(input_mod_3 == "Iwamoto_WDD3") or (input_mod_3 == "Iwamoto_CDD1") or \
			(input_mod_3 == "Iwamoto_CDD2") or \
			(input_mod_3 == "Badenes_DDTa") or (input_mod_3 == "Badenes_DDTb") or \
			(input_mod_3 == "Badenes_DDTc") or (input_mod_3 == "Badenes_DDTd") or \
			(input_mod_3 == "Badenes_DDTe") or (input_mod_3 == "Badenes_DDTf") or \
			(input_mod_3 == "Badenes_DDTa_undecayed") or \
			(input_mod_3 == "Badenes_DDTb_undecayed") or \
			(input_mod_3 == "Badenes_DDTc_undecayed") or \
			(input_mod_3 == "Badenes_DDTd_undecayed") or \
			(input_mod_3 == "Badenes_DDTe_undecayed") or \
			(input_mod_3 == "Badenes_DDTf_undecayed") or \
			(input_mod_3 == "Maeda_C-DEF") or \
			(input_mod_3 == "Maeda_C-DDT") or \
			(input_mod_3 == "Maeda_O-DDT") or \
			(input_mod_3 == "Waldman_CO.45HE.2") or \
			(input_mod_3 == "Waldman_CO.55HE.2") or \
			(input_mod_3 == "Waldman_CO.5HE.15") or \
			(input_mod_3 == "Waldman_CO.5HE.2") or \
			(input_mod_3 == "Waldman_CO.5HE.2C.3") or \
			(input_mod_3 == "Waldman_CO.5HE.2N.02") or \
			(input_mod_3 == "Waldman_CO.5HE.3") or \
			(input_mod_3 == "Waldman_CO.6HE.2"):
			model3 = sum_SNIa(input_mod_3,elements)
		elif (input_mod_3 == "0") or (input_mod_3 == "0.001") or \
			(input_mod_3 == "0.004") or (input_mod_3 == "0.02") or \
			(input_mod_3 == "0_undecayed") or (input_mod_3 == "0.001_undecayed") or \
			(input_mod_3 == "0.004_undecayed") or (input_mod_3 == "0.02_undecayed"):
			model3 = integrate_SNcc(input_mod_3,alpha,elem_name)
		else:
			print "Error: The model", input_mod_3, "does not exit, or could not be found."
			sys.exit()


	if number_of_models>3:
		if (input_mod_4 == "Iwamoto_W7") or (input_mod_4 == "Iwamoto_W70") or \
			(input_mod_4 == "Iwamoto_WDD1") or (input_mod_4 == "Iwamoto_WDD2") or \
			(input_mod_4 == "Iwamoto_WDD3") or (input_mod_4 == "Iwamoto_CDD1") or \
			(input_mod_4 == "Iwamoto_CDD2") or \
			(input_mod_4 == "Badenes_DDTa") or (input_mod_4 == "Badenes_DDTb") or \
			(input_mod_4 == "Badenes_DDTc") or (input_mod_4 == "Badenes_DDTd") or \
			(input_mod_4 == "Badenes_DDTe") or (input_mod_4 == "Badenes_DDTf") or \
			(input_mod_4 == "Badenes_DDTa_undecayed") or \
			(input_mod_4 == "Badenes_DDTb_undecayed") or \
			(input_mod_4 == "Badenes_DDTc_undecayed") or \
			(input_mod_4 == "Badenes_DDTd_undecayed") or \
			(input_mod_4 == "Badenes_DDTe_undecayed") or \
			(input_mod_4 == "Badenes_DDTf_undecayed") or \
			(input_mod_4 == "Maeda_C-DEF") or \
			(input_mod_4 == "Maeda_C-DDT") or \
			(input_mod_4 == "Maeda_O-DDT") or \
			(input_mod_4 == "Waldman_CO.45HE.2") or \
			(input_mod_4 == "Waldman_CO.55HE.2") or \
			(input_mod_4 == "Waldman_CO.5HE.15") or \
			(input_mod_4 == "Waldman_CO.5HE.2") or \
			(input_mod_4 == "Waldman_CO.5HE.2C.3") or \
			(input_mod_4 == "Waldman_CO.5HE.2N.02") or \
			(input_mod_4 == "Waldman_CO.5HE.3") or \
			(input_mod_4 == "Waldman_CO.6HE.2"):
			model4 = sum_SNIa(input_mod_4,elements)
		elif (input_mod_4 == "0") or (input_mod_4 == "0.001") or \
			(input_mod_4 == "0.004") or (input_mod_4 == "0.02") or \
			(input_mod_4 == "0_undecayed") or (input_mod_4 == "0.001_undecayed") or \
			(input_mod_4 == "0.004_undecayed") or (input_mod_4 == "0.02_undecayed"):
			model4 = integrate_SNcc(input_mod_4,alpha,elem_name)
		else:
			print "Error: The model", input_mod_4, "does not exit, or could not be found."
			sys.exit()




	print "Model: ", input_mod_1
	print model1
	if number_of_models>1:
		print "Model: ", input_mod_2
		print model2
	if number_of_models>2:
		print "Model: ", input_mod_3
		print model3
	if number_of_models>3:
		print "Model: ", input_mod_4
		print model4


	x=numpy.arange(len(elements))
	x=x*0.0
	if number_of_models>1:
		y=numpy.arange(len(elements))
		y=y*0.0
	if number_of_models>2:
		z=numpy.arange(len(elements))
		z=z*0.0
	if number_of_models>3:
		w=numpy.arange(len(elements))
		w=w*0.0

	for i in range(len(elements)):
		x[i]=model1[i]/(mass[i]*lodders[i])
		if number_of_models>1:
			y[i]=model2[i]/(mass[i]*lodders[i])
		if number_of_models>2:
			z[i]=model3[i]/(mass[i]*lodders[i])
		if number_of_models>3:
			w[i]=model4[i]/(mass[i]*lodders[i])

	




	# Fit...
	SNe=numpy.arange(number_of_models)
	x0=SNe

	if number_of_models==1:
		modmatrix = numpy.vstack((x))
		#print optimization.leastsq(fitfunc1, x0, args=(modmatrix, data))
		bestSNe=numpy.array([0])
		bestSNe[0]=fitfunc1(modmatrix, data)
		bestfit=numpy.dot(modmatrix,bestSNe)

	if number_of_models==2:
		modmatrix = numpy.vstack((x,y))
		#print optimization.leastsq(fitfunc2, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc2, x0, args=(modmatrix, data))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)

	if number_of_models==3:
		modmatrix = numpy.vstack((x,y,z))
		#print optimization.leastsq(fitfunc3, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc3, x0, args=(modmatrix, data))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)

	if number_of_models==4:
		modmatrix = numpy.vstack((x,y,z,w))
		#print optimization.leastsq(fitfunc4, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc4, x0, args=(modmatrix, data))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)






	#Compute chi2...
	chi = 0.0
	print "Optimal parameters:", bestSNe
	for i in range(len(elements)):
		if number_of_models==1:
			chi = chi + (data[i] - bestSNe[0]*x[i])**2 / errors[i]**2
		if number_of_models==2:
			chi = chi + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i])**2 / errors[i]**2
		if number_of_models==3:
			chi = chi + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i])**2 / errors[i]**2
		if number_of_models==4:
			chi = chi + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i] - bestSNe[3]*w[i])**2 / errors[i]**2

	print "chi^2 / d.o.f. = ", float(chi),'/',(len(elements)-number_of_models), "=", chi / (len(elements)-number_of_models)

	chiout = open("currentchi.txt", 'w+')
	chiout.write("%s" % int(chi*100000))		# Used for serial_fitter.sh. *100000-> to distinguish the decimals (will be corrected in sort_summary_models.py)





	#Display QDP to screen...
	print "--------------------------QDP FILE--------------------------" 
	print "skip single"
	print "line on 2,3,4,5,6"
	print "mark 0 on 1"
	print "r y -0.1 3"
	print "la x Atomic Number"
	print "la y Abundance (proto-solar)"
	print "READ Serr 1,2"
	print "! Abundance measurements"
	for el in range(len(elements)):
		print elements[el], 0, data[el], errors[el]
	print "NO"
	print "! Total of models"
	for el in range(len(elements)):
		print elements[el], 0, bestfit[el], 0
	print "NO"
	print "! Model:", input_mod_1
	for el in range(len(elements)):
		print elements[el], 0, bestSNe[0]*x[el], 0
	if number_of_models>1:
		print "NO"
		print "! Model:", input_mod_2
		for el in range(len(elements)):
			print elements[el], 0, bestSNe[1]*y[el], 0
	if number_of_models>2:
		print "NO"
		print "! Model:", input_mod_3
		for el in range(len(elements)):
			print elements[el], 0, bestSNe[2]*z[el], 0
	if number_of_models>3:
		print "NO"
		print "! Model:", input_mod_4
		for el in range(len(elements)):
			print elements[el], 0, bestSNe[3]*w[el], 0





	#Write into QDP file...
	qdpout = open("output.qdp", 'w+')
	qdpout.write("skip single\n")
	qdpout.write("line on 2,3,4,5,6\n")
	qdpout.write("mark 0 on 1\n")
	qdpout.write("r y -0.1 3\n")
	qdpout.write("READ Serr 1,2\n")
	qdpout.write("! Abundance measurements\n")

	for el in range(len(elements)):
		qdpout.write("%s %s %s %s\n" % (elements[el], 0, data[el], errors[el]))

	qdpout.write("NO\n")
	qdpout.write("! Total of models\n")

	for el in range(len(elements)):
		qdpout.write("%s %s %s %s\n" % (elements[el], 0, bestfit[el], 0))

	qdpout.write("NO\n")
	qdpout.write("! Model: %s\n" % input_mod_1)

	for el in range(len(elements)):
		qdpout.write("%s %s %s %s\n" % (elements[el], 0, bestSNe[0]*x[el], 0))

	if number_of_models>1:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_2)

		for el in range(len(elements)):
			qdpout.write("%s %s %s %s\n" % (elements[el], 0, bestSNe[1]*y[el], 0))

	if number_of_models>2:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_3)

		for el in range(len(elements)):
			qdpout.write("%s %s %s %s\n" % (elements[el], 0, bestSNe[2]*z[el], 0))

	if number_of_models>3:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_4)

		for el in range(len(elements)):
			qdpout.write("%s %s %s %s\n" % (elements[el], 0, bestSNe[3]*w[el], 0))








	#Write into "barplot.dat"...
	barplout = open("barplot.dat", 'w+')

	if number_of_models==1:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(elements[el]), data[el], errors[el], bestSNe[0]*x[el], 0.0, 0.0))

	if number_of_models==2:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(elements[el]), data[el], errors[el], bestSNe[0]*x[el], bestSNe[1]*y[el], 0.0))

	if number_of_models==3:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(elements[el]), data[el], errors[el], bestSNe[0]*x[el], bestSNe[1]*y[el], bestSNe[2]*z[el]))

	if number_of_models==4:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s %s\n" % (int(elements[el]), data[el], errors[el], bestSNe[0]*x[el], bestSNe[1]*y[el], bestSNe[2]*z[el], bestSNe[3]*w[el]))






















# # # # # # # FUNCTIONS # # # # # # #
# # # # # # # # # # # # # # # # # # #



## Integrates the SNcc products over a mass range of 13-40 M_sun
## and assuming a certain IMF.
def integrate_SNcc( metal, alpha , elem_name ):

	m=[13,15,18,20,25,30,40]
	dm=[2,2.5,2.5,3.5,5,7.5,10]
	y=numpy.arange(len(elem_name))
	y=y*0.0

	alpha = float(alpha)
	#alpha=-2.35 (Salpeter)
	#alpha=-1.0 (Top)
    
  #for z in metal:     
	n=0

	for el in elem_name:
		data_mod=numpy.loadtxt('SNcc/nomoto/'+str(metal)+'/'+str(el)+'.txt') #Open file
		yiel=numpy.arange(len(m))
		yiel=0.0*yiel
		for j in range(0,len(data_mod)):       #Merge the isotopes of a same element
			yiel = yiel + data_mod[j,]

		sumt=0.
  		sumn=0.
	    
  		for i in numpy.arange(len(m)):    #Integration (see de Plaa+ 07)
			sumt=sumt+yiel[i]*dm[i]*m[i]**(alpha)
			sumn=sumn+dm[i]*m[i]**(alpha)


  		y[n]=sumt/sumn
  		n=n+1

	return y





## Sums the SNIa products.
def sum_SNIa( SNIa_model , elements ):

	y=numpy.arange(len(elements))
	y=y*0.0

	data_mod=numpy.loadtxt('SNIa/'+str(SNIa_model)+'.txt')

	for el in range(len(elements)):
		for i in range(len(data_mod)):
			if data_mod[i,0] == float(elements[el]):
				y[el] = y[el] + data_mod[i,2]

	return y















## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).
def fitfunc1(modmatrix, data):
	Fe_index=-1
	for i in range(len(data)):	# Calculate the array index for Fe
		if data[i] ==1:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: You must set the Fe abundance to 1 in your input file. (The models are calculated relatively to Fe).")

	SNe=1.0/modmatrix[Fe_index,0]
	return SNe



## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).
def fitfunc2(SNe,modmatrix, data):
	Fe_index=-1
	for i in range(len(data)):	# Calculate the array index for Fe
		if data[i] ==1:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: You must set the Fe abundance to 1 in your input file. (The models are calculated relatively to Fe).")

	SNe[1]=(data[Fe_index]-SNe[0]*modmatrix[0,Fe_index])/modmatrix[1,Fe_index]
	return (data - numpy.dot(SNe,modmatrix))




## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).

def fitfunc3(SNe,modmatrix, data):
	Fe_index=-1
	for i in range(len(data)):	# Calculate the array index for Fe
		if data[i] ==1:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: You must set the Fe abundance to 1 in your input file. (The models are calculated relatively to Fe).")

	SNe[2]=(data[Fe_index]-SNe[0]*modmatrix[0,Fe_index]-SNe[1]*modmatrix[1,Fe_index])/modmatrix[2,Fe_index]
	return (data - numpy.dot(SNe,modmatrix))





## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).

def fitfunc4(SNe,modmatrix, data):
	Fe_index=-1
	for i in range(len(data)):	# Calculate the array index for Fe
		if data[i] ==1:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: You must set the Fe abundance to 1 in your input file. (The models are calculated relatively to Fe).")

	SNe[3]=(data[Fe_index]-SNe[0]*modmatrix[0,Fe_index]-SNe[1]*modmatrix[1,Fe_index]-SNe[2]*modmatrix[2,Fe_index])/modmatrix[3,Fe_index]
	return (data - numpy.dot(SNe,modmatrix))


















#============= Gamma, Incomplete Gamma ===========

def gammln(xx):
    #"""Logarithm of the gamma function."""
    global gammln_cof, gammln_stp
    rt2 = math.sqrt(2.)
    gammln_cof = numpy.array([76.18009173, -86.50532033, 24.01409822, -1.231739516e0, 0.120858003e-2, -0.536382e-5])
    gammln_stp = 2.50662827465
    x = xx - 1.
    tmp = x + 5.5
    tmp = (x + 0.5)*math.log(tmp) - tmp
    ser = 1.
    for j in range(6):
        x = x + 1.
        ser = ser + gammln_cof[j]/x
    return tmp + math.log(gammln_stp*ser)

def gser(a, x, itmax=700, eps=3.e-7):
    #"""Series approx'n to the incomplete gamma function."""
    gln = gammln(a)
    if (x < 0.):
        raise bad_arg, x
    if (x == 0.):
        return(0.)
    ap = a
    sum = 1. / a
    delta = sum
    n = 1
    while n <= itmax:
        ap = ap + 1.
        delta = delta * x / ap
        sum = sum + delta
        if (abs(delta) < abs(sum)*eps):
            return (sum * math.exp(-x + a*math.log(x) - gln), gln)
        n = n + 1
    raise max_iters, str((abs(delta), abs(sum)*eps))


def gcf(a, x, itmax=200, eps=3.e-7):
    #"""Continued fraction approx'n of the incomplete gamma function."""
    gln = gammln(a)
    gold = 0.
    a0 = 1.
    a1 = x
    b0 = 0.
    b1 = 1.
    fac = 1.
    n = 1
    while n <= itmax:
        an = n
        ana = an - a
        a0 = (a1 + a0*ana)*fac
        b0 = (b1 + b0*ana)*fac
        anf = an*fac
        a1 = x*a0 + anf*a1
        b1 = x*b0 + anf*b1
        if (a1 != 0.):
            fac = 1. / a1
            g = b1*fac
            if (abs((g-gold)/g) < eps):
                return (g*math.exp(-x+a*math.log(x)-gln), gln)
            gold = g
        n = n + 1
    raise max_iters, str(abs((g-gold)/g))


def gammp(a, x):
    #"""Incomplete gamma function."""
    if (x < 0. or a <= 0.):
        raise ValueError, (a, x)
    if (x < a+1.):
        return gser(a,x)[0]
    else:
        return 1.-gcf(a,x)[0]

def gammq(a, x):
    #"""Incomplete gamma function."""
    if (x < 0. or a <= 0.):
        raise ValueError, repr((a, x))
    if (x < a+1.):
        return 1.-gser(a,x)[0]
    else:
        return gcf(a,x)[0]

def factrl(n, ntop=0, prev=numpy.ones((33),dtype=float)):
#"""Factorial of n. The first 33 values are stored as they are calculated to speed up subsequent calculations."""
    if n < 0:
        raise ValueError, 'Negative argument!'
    elif n <= ntop:
        return prev[n]
    elif n <= 32:
        for j in range(ntop+1, n+1):
            prev[j] = j * prev[j-1]
            ntop = n
        return prev[n]
    else:
        return math.exp(gammln(n+1.))

def factln(n, prev=numpy.array(101*(-1.,))):
#"""Log factorial of n. Values for n=0 to 100 are stored as they are calculated to speed up subsequent calls."""
    if n < 0:
        raise ValueError, 'Negative argument!'
    elif n <= 100:
        if prev[n] < 0:
           prev[n] = gammln(n+1.)
        return prev[n]
    else:
        return gammln(n+1.)

def combln(Ntot, n):
    #"""Log of number of combinations of n items from Ntot total."""
    return factln(Ntot) - factln(n) - factln(Ntot-n)

   
main()
