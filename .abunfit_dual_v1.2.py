#!/usr/bin/python

import os
import sys
import math
import numpy
import scipy.optimize as optimization


def main():

	if len(sys.argv) < 4:
		sys.exit("Error: Uncorrect arguments.\n\n \
		PROGRAM: abunfit_dual.py\n\n \
		SYNOPSIS: ./abunfit_dual.py [input_file] [number_of_models] [alpha] [SNe_model_1] [SNe_model_2] ... \n\n \
		DESCRIPTION: Fits a set of intra-cluster abundances with a combination of \n \
		supernovae (SNe) yields models. Several models (either SNcc or SNIa) are available \n \
		and can be updated or created. The models are fitted simultaneously with proto-solar values\n \
		(and their uncertainties), in order to account for *both* the galactic and the ICM enrichments.\n\n \
		ARGUMENTS: \n\n \
		        [input_file] -- (Location of the) Text input file of the set of measured  \n \
		                        ICM abundances. This file must contain 3 columns:  \n \
		                        (1) the element number (Z), \n \
		                        (2) the measured abundance, and \n \
		                        (3) its associated (symmetrical) uncertainties.\n\n \
		        [number_of_models] -- Integer (from 1 to 4) indicating the number of \n \
 			                      models to be fitted simultaneously (min: 2, max: 4).\n \
		        [alpha] -- Slope index of the assumed initial mass function (IMF). \n \
		                   Ex: for Salpeter IMF, alpha = -2.35.\n \
		        [model<1/2/...>] -- Name of the model. Currently available models: \n \
		                            ==============\n \
		                            ==== SNIa ====\n \
		                            ==============\n \
		                            == Iwamoto et al. (1999)\n \
		                             - Iwamoto_W7\n \
		                             - Iwamoto_W70\n \
		                             - Iwamoto_WDD1\n \
		                             - Iwamoto_WDD2\n \
		                             - Iwamoto_WDD3\n \
		                             - Iwamoto_CDD1\n \
		                             - Iwamoto_CDD2\n \
		                            == Badenes et al. (2006)\n \
		                               (Based on Tycho's best spectral fits...)\n \
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
		                            == Seitenzahl et al. (2013)\n \
		                               (3-D delayed-detonation explosion models...)\n \
		                             - Seitenzahl_N1\n \
		                             - Seitenzahl_N3\n \
		                             - Seitenzahl_N5\n \
		                             - Seitenzahl_N10\n \
		                             - Seitenzahl_N20\n \
		                             - Seitenzahl_N40\n \
		                             - Seitenzahl_N100H\n \
		                             - Seitenzahl_N100\n \
		                             - Seitenzahl_N100L\n \
		                             - Seitenzahl_N150\n \
		                             - Seitenzahl_N200\n \
		                             - Seitenzahl_N300C\n \
		                             - Seitenzahl_N1600\n \
		                             - Seitenzahl_N1600C\n \
		                             - Seitenzahl_N100_Z0.5\n \
		                             - Seitenzahl_N100_Z0.1\n \
		                             - Seitenzahl_N100_Z0.01\n \
		                            == Fink et al. (2014)\n \
		                               (3-D deflagration explosion models...)\n \
		                             - Fink_N1def\n \
		                             - Fink_N3def\n \
		                             - Fink_N5def\n \
		                             - Fink_N10def\n \
		                             - Fink_N20def\n \
		                             - Fink_N40def\n \
		                             - Fink_N100Hdef\n \
		                             - Fink_N100def\n \
		                             - Fink_N100Ldef\n \
		                             - Fink_N150def\n \
		                             - Fink_N200def\n \
		                             - Fink_N300Cdef\n \
		                             - Fink_N1600def\n \
		                             - Fink_N1600Cdef\n \
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
		                            ==============\n \
		                            ==== SNcc ====\n \
		                            ==============\n \
		                            == Chieffi & Limongi (2004)\n \
		                               (One initial metallicity per model...)\n \
		                             - Chieffi_0\n \
		                             - Chieffi_1E-6\n \
		                             - Chieffi_1E-4\n \
		                             - Chieffi_1E-3\n \
		                             - Chieffi_6E-3\n \
		                             - Chieffi_2E-2\n \
		                            == Nomoto et al. (2006)\n \
		                               (One initial metallicity per model...)\n \
		                             - Nomoto_0\n \
		                             - Nomoto_0.001\n \
		                             - Nomoto_0.004\n \
		                             - Nomoto_0.02\n \
		                            == Nomoto et al. (2013)\n \
		                               (One initial metallicity per model...)\n \
		                             = SNcc: 11 - 40 (100) M_sun\n \
		                              - N13_SNcc_0\n \
		                              - N13_SNcc_0.001\n \
		                              - N13_SNcc_0.004\n \
		                              - N13_SNcc_0.008\n \
		                              - N13_SNcc_0.02\n \
		                              - N13_SNcc_0.05\n \
		                             = PISNe: 140 - 300 M_sun\n \
		                              - N13_PISNe_0\n \
		                             = All SNe (SNcc+PISNe): 11 - 300 M_sun\n \
		                              - N13_SNe_0\n\n \
		EXAMPLE: ./abunfit_dual.py cluster/average_abun_clusters.out 2 -2.35 Iwamoto_WDD2 Nomoto_0.02\n\n \
		VERSION: 1.2 (February 2016)\n \
		         1.1 (January 2016)\n \
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
	alpha = float(alpha)


	# Read proto-solar data.
	f4 = file("proto-solar.out")
	sol_abun = numpy.loadtxt(f4)

	sol_elements = sol_abun[:,0]
	sol_data = sol_abun[:,1]
	sol_errors = sol_abun[:,2]


	# Give to each considered element its respective name.
	cl_elem_name = elem_name(elements)
	sol_elem_name = elem_name(sol_elements)






	# Read atomic masses. 
	# mass[i] = atomic mass of the ith element.
	f2 = file("atomic_masses.txt")
	atom_mass = numpy.loadtxt(f2)

	mass = numpy.arange(len(elements))
	mass = mass*0.0
	for i in range(len(elements)):
		mass[i] = atom_mass[elements[i]-1,1]

	sol_mass = numpy.arange(len(sol_elements))		### Proto-solar fits
	sol_mass = sol_mass*0.0					### Proto-solar fits
	for i in range(len(sol_elements)):			### Proto-solar fits
		sol_mass[i] = atom_mass[sol_elements[i]-1,1]	### Proto-solar fits





	# Read abundance tables (Lodders et al. 2009). 
	# lodders[i] = protosolar abun. of the ith element (rel. to Si==10^6).
	f3 = file("lodders09.txt")
	abun_table = numpy.loadtxt(f3)

	lodders = numpy.arange(len(elements))
	lodders = lodders*0.0
	for i in range(len(elements)):
		lodders[i] = abun_table[elements[i]-1,1]

	sol_lodders = numpy.arange(len(sol_elements))			### Proto-solar fits
	sol_lodders = sol_lodders*0.0					### Proto-solar fits
	for i in range(len(sol_elements)):				### Proto-solar fits
		sol_lodders[i] = abun_table[sol_elements[i]-1,1]	### Proto-solar fits






	# Read list of available SNe models. 
	f4 = file("list_models_SNIa.txt")
	lines = f4.readlines()
	list_models_SNIa = [remove_end_char(s) for s in lines]
	f5 = file("list_models_SNcc.txt")
	lines = f5.readlines()
	list_models_SNcc = [remove_end_char(s) for s in lines]






	# Load SNIa and/or SNcc models...

	if (input_mod_1 in list_models_SNIa):
		model1 = sum_SNIa(input_mod_1,elements)
		sol_model1 = sum_SNIa(input_mod_1,sol_elements)
	elif (input_mod_1 in list_models_SNcc):
		model1 = integrate_SNcc(input_mod_1,alpha,cl_elem_name)
		sol_model1 = integrate_SNcc(input_mod_1,alpha,sol_elem_name)
	else:
		print "Error: The model", input_mod_1, "does not exist, or could not be found."
		sys.exit()


	if number_of_models>1:
		if (input_mod_2 in list_models_SNIa):
			model2 = sum_SNIa(input_mod_2,elements)
			sol_model2 = sum_SNIa(input_mod_2,sol_elements)
		elif (input_mod_2 in list_models_SNcc):
			model2 = integrate_SNcc(input_mod_2,alpha,cl_elem_name)
			sol_model2 = integrate_SNcc(input_mod_2,alpha,sol_elem_name)
		else:
			print "Error: The model", input_mod_2, "does not exist, or could not be found."
			sys.exit()


	if number_of_models>2:
		if (input_mod_3 in list_models_SNIa):
			model3 = sum_SNIa(input_mod_3,elements)
			sol_model3 = sum_SNIa(input_mod_3,sol_elements)
		elif (input_mod_3 in list_models_SNcc):
			model3 = integrate_SNcc(input_mod_3,alpha,cl_elem_name)
			sol_model3 = integrate_SNcc(input_mod_3,alpha,sol_elem_name)
		else:
			print "Error: The model", input_mod_3, "does not exist, or could not be found."
			sys.exit()


	if number_of_models>3:
		if (input_mod_4 in list_models_SNIa):
			model4 = sum_SNIa(input_mod_4,elements)
			sol_model4 = sum_SNIa(input_mod_4,sol_elements)
		elif (input_mod_4 in list_models_SNcc):
			model4 = integrate_SNcc(input_mod_4,alpha,cl_elem_name)
			sol_model4 = integrate_SNcc(input_mod_4,alpha,sol_elem_name)
		else:
			print "Error: The model", input_mod_4, "does not exist, or could not be found."
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
	sol_x=numpy.arange(len(sol_elements))
	sol_x=sol_x*0.0
	if number_of_models>1:
		y=numpy.arange(len(elements))
		y=y*0.0
		sol_y=numpy.arange(len(sol_elements))
		sol_y=sol_y*0.0
	if number_of_models>2:
		z=numpy.arange(len(elements))
		z=z*0.0
		sol_z=numpy.arange(len(sol_elements))
		sol_z=sol_z*0.0
	if number_of_models>3:
		w=numpy.arange(len(elements))
		w=w*0.0
		sol_w=numpy.arange(len(sol_elements))
		sol_w=sol_w*0.0


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=model1[i]/(mass[i]*lodders[i])
		if number_of_models>1:
			y[i]=model2[i]/(mass[i]*lodders[i])
		if number_of_models>2:
			z[i]=model3[i]/(mass[i]*lodders[i])
		if number_of_models>3:
			w[i]=model4[i]/(mass[i]*lodders[i])

	# Define X, Y (, Z, W): abun / supernova	###Proto-solar fit
	for i in range(len(sol_elements)):
		sol_x[i]=sol_model1[i]/(sol_mass[i]*sol_lodders[i])
		if number_of_models>1:
			sol_y[i]=sol_model2[i]/(sol_mass[i]*sol_lodders[i])
		if number_of_models>2:
			sol_z[i]=sol_model3[i]/(sol_mass[i]*sol_lodders[i])
		if number_of_models>3:
			sol_w[i]=sol_model4[i]/(sol_mass[i]*sol_lodders[i])

	




	# Fit...
	SNe=numpy.arange(number_of_models)
	x0=SNe
	sol_SNe=numpy.arange(number_of_models)	###Proto-solar fit
	sol_x0=sol_SNe				###Proto-solar fit

	if number_of_models==1:
		modmatrix = numpy.vstack((x))
		sol_modmatrix = numpy.vstack((sol_x))
		#print optimization.leastsq(fitfunc1, x0, args=(modmatrix, data))
		bestSNe=numpy.array([0])
		bestSNe[0]=fitfunc1(modmatrix, cl_abun)
		bestfit=numpy.dot(modmatrix,bestSNe)
		sol_bestSNe=numpy.array([0])				###Proto-solar fit
		sol_bestSNe[0]=fitfunc1(sol_modmatrix, sol_abun)	###Proto-solar fit
		sol_bestfit=numpy.dot(sol_modmatrix,sol_bestSNe)	###Proto-solar fit

	if number_of_models==2:
		modmatrix = numpy.vstack((x,y))
		sol_modmatrix = numpy.vstack((sol_x,sol_y))
		#print optimization.leastsq(fitfunc2, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc2, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)
		sol_bestSNe=optimization.leastsq(fitfunc2, x0, args=(sol_modmatrix, sol_abun))[0]	###Proto-solar fit
		sol_bestfit=numpy.dot(sol_bestSNe,sol_modmatrix)					###Proto-solar fit

	if number_of_models==3:
		modmatrix = numpy.vstack((x,y,z))
		sol_modmatrix = numpy.vstack((sol_x,sol_y,sol_z))
		#print optimization.leastsq(fitfunc3, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc3, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)
		sol_bestSNe=optimization.leastsq(fitfunc3, x0, args=(sol_modmatrix, sol_abun))[0]	###Proto-solar fit
		sol_bestfit=numpy.dot(sol_bestSNe,sol_modmatrix)					###Proto-solar fit

	if number_of_models==4:
		modmatrix = numpy.vstack((x,y,z,w))
		sol_modmatrix = numpy.vstack((sol_x,sol_y,sol_z,sol_w))
		#print optimization.leastsq(fitfunc4, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc4, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)
		sol_bestSNe=optimization.leastsq(fitfunc4, x0, args=(sol_modmatrix, sol_abun))[0]	###Proto-solar fit
		sol_bestfit=numpy.dot(sol_bestSNe,sol_modmatrix)					###Proto-solar fit








	for i in range(len(bestSNe)):
		if bestSNe[i] < 0.0:
			print "Error: One model has a negative contribution to the enrichment. This is physically not allowed."
			chiout = open("currentchi.txt", 'w+')
			chiout.write("%s" % int(99999999999))		# Used for serial_fitter.sh. *100000-> to distinguish the decimals (will be corrected in sort_summary_models.py)
			sys.exit()

	for i in range(len(sol_bestSNe)):
		if sol_bestSNe[i] < 0.0:
			print "Error: One model has a negative contribution to the enrichment. This is physically not allowed."
			chiout = open("currentchi.txt", 'w+')
			chiout.write("%s" % int(99999999999))		# Used for serial_fitter.sh. *100000-> to distinguish the decimals (will be corrected in sort_summary_models.py)
			sys.exit()



	#Compute chi2...
	chi = 0.0
	chi_a = 0.0
	chi_b = 0.0
	print "Optimal parameters (ICM):", bestSNe
	print "Optimal parameters (proto-solar):", sol_bestSNe
	for i in range(len(elements)):
		if number_of_models==1:
			chi_a = chi_a + (data[i] - bestSNe[0]*x[i])**2 / errors[i]**2
		if number_of_models==2:
			chi_a = chi_a + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i])**2 / errors[i]**2
		if number_of_models==3:
			chi_a = chi_a + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i])**2 / errors[i]**2
		if number_of_models==4:
			chi_a = chi_a + (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i] - bestSNe[3]*w[i])**2 / errors[i]**2

	for i in range(len(sol_elements)):
		if number_of_models==1:
			chi_b = chi_b + (sol_data[i] - sol_bestSNe[0]*sol_x[i])**2 / sol_errors[i]**2
		if number_of_models==2:
			chi_b = chi_b + (sol_data[i] - sol_bestSNe[0]*sol_x[i] - sol_bestSNe[1]*sol_y[i])**2 / sol_errors[i]**2
		if number_of_models==3:
			chi_b = chi_b + (sol_data[i] - sol_bestSNe[0]*sol_x[i] - sol_bestSNe[1]*sol_y[i] - sol_bestSNe[2]*sol_z[i])**2 / sol_errors[i]**2
		if number_of_models==4:
			chi_b = chi_b + (sol_data[i] - sol_bestSNe[0]*sol_x[i] - sol_bestSNe[1]*sol_y[i] - sol_bestSNe[2]*sol_z[i] - sol_bestSNe[3]*sol_w[i])**2 / sol_errors[i]**2

	chi = chi_a + chi_b

	print "chi^2 / d.o.f. = ", float(chi),'/',(len(elements) + len(sol_elements)- 2*number_of_models), "=", \
	       chi / (len(elements) + len(sol_elements)- 2*number_of_models)

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





	#Write into QDP file... ### Proto-solar fit
	qdpout = open("output_sol.qdp", 'w+')
	qdpout.write("skip single\n")
	qdpout.write("line on 2,3,4,5,6\n")
	qdpout.write("mark 0 on 1\n")
	qdpout.write("r y -0.1 3\n")
	qdpout.write("READ Serr 1,2\n")
	qdpout.write("! Abundance measurements\n")

	for el in range(len(sol_elements)):
		qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_data[el], sol_errors[el]))

	qdpout.write("NO\n")
	qdpout.write("! Total of models\n")

	for el in range(len(sol_elements)):
		qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_bestfit[el], 0))

	qdpout.write("NO\n")
	qdpout.write("! Model: %s\n" % input_mod_1)

	for el in range(len(sol_elements)):
		qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_bestSNe[0]*sol_x[el], 0))

	if number_of_models>1:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_2)

		for el in range(len(sol_elements)):
			qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_bestSNe[1]*sol_y[el], 0))

	if number_of_models>2:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_3)

		for el in range(len(sol_elements)):
			qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_bestSNe[2]*sol_z[el], 0))

	if number_of_models>3:
		qdpout.write("NO\n")
		qdpout.write("! Model: %s\n" % input_mod_4)

		for el in range(len(sol_elements)):
			qdpout.write("%s %s %s %s\n" % (sol_elements[el], 0, sol_bestSNe[3]*sol_w[el], 0))








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




	#Write into "barplot.dat"... ### Proto-solar fit
	barplout = open("barplot_sol.dat", 'w+')

	if number_of_models==1:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(sol_elements[el]), sol_data[el], sol_errors[el], sol_bestSNe[0]*x[el], 0.0, 0.0))

	if number_of_models==2:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(sol_elements[el]), sol_data[el], sol_errors[el], sol_bestSNe[0]*x[el], sol_bestSNe[1]*y[el], 0.0))

	if number_of_models==3:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s\n" % (int(sol_elements[el]), sol_data[el], sol_errors[el], sol_bestSNe[0]*x[el], sol_bestSNe[1]*y[el], sol_bestSNe[2]*z[el]))

	if number_of_models==4:
		for el in range(len(elements)):
			barplout.write("%s %s %s %s %s %s %s\n" % (int(sol_elements[el]), sol_data[el], sol_errors[el], sol_bestSNe[0]*x[el], sol_bestSNe[1]*y[el], sol_bestSNe[2]*z[el], sol_bestSNe[3]*w[el]))



	#Last check, to be sure that Fe is normalized to 1.0 ...
	checkFe(cl_abun)

















# # # # # # # FUNCTIONS # # # # # # #
# # # # # # # # # # # # # # # # # # #



## Integrates the SNcc products over a mass range of 13-40 M_sun
## and assuming a certain IMF.
def integrate_SNcc( input_model, alpha , elem_name ):

	if (input_model == "Nomoto_0") or (input_model == "Nomoto_0.001") or \
		(input_model == "Nomoto_0.004") or (input_model == "Nomoto_0.02") or \
		(input_model == "Nomoto_0_undecayed") or (input_model == "Nomoto_0.001_undecayed") or \
		(input_model == "Nomoto_0.004_undecayed") or (input_model == "Nomoto_0.02_undecayed"):
		m=[13,15,18,20,25,30,40]
		dm=[2,2.5,2.5,3.5,5,7.5,10]
	elif (input_model == "Chieffi_0") or (input_model == "Chieffi_1E-6") or \
		(input_model == "Chieffi_1E-4") or (input_model == "Chieffi_1E-3") or \
		(input_model == "Chieffi_6E-3") or (input_model == "Chieffi_2E-2"):
		m=[13,15,20,25,30,35]
		dm=[2,3.5,5,5,5,10]
	elif (input_model == "N13_SNcc_0"):
		m=[11,13,15,18,20,25,30,40,100,140]
		dm=[1.5,2,2.5,2.5,3.5,5,7.5,35,50,20]
	elif (input_model == "N13_SNcc_0.001") or (input_model == "N13_SNcc_0.004") or \
		(input_model == "N13_SNcc_0.008") or (input_model == "N13_SNcc_0.02") or \
		(input_model == "N13_SNcc_0.05"):
		m=[13,15,18,20,25,30,40]
		dm=[2,2.5,2.5,3.5,5,7.5,10]
	elif (input_model == "N13_PISNe_0"):
		m=[140,150,170,200,270,300]
		dm=[5,15,25,50,50,15]
	elif (input_model == "N13_SNe_0"):
		m=[11,13,15,18,20,25,30,40,100,140,150,170,200,270,300]
		dm=[1.5,2,2.5,2.5,3.5,5,7.5,35,50,25,15,25,50,50,15]
	y=numpy.arange(len(elem_name))
	y=y*0.0

	alpha = float(alpha)
	#alpha=-2.35 (Salpeter)
	#alpha=-1.0 (Top)
    
  #for z in metal:     
	n=0

	for el in elem_name:
		if (input_model == "Nomoto_0"):
			data_mod=numpy.loadtxt('SNcc/nomoto/0/'+str(el)+'.txt') #Open file
		elif (input_model == "Nomoto_0.001"):
			data_mod=numpy.loadtxt('SNcc/nomoto/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "Nomoto_0.004"):
			data_mod=numpy.loadtxt('SNcc/nomoto/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "Nomoto_0.02"):
			data_mod=numpy.loadtxt('SNcc/nomoto/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_0"):
			data_mod=numpy.loadtxt('SNcc/chieffi/0/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_1E-6"):
			data_mod=numpy.loadtxt('SNcc/chieffi/1E-6/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_1E-4"):
			data_mod=numpy.loadtxt('SNcc/chieffi/1E-4/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_1E-3"):
			data_mod=numpy.loadtxt('SNcc/chieffi/1E-3/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_6E-3"):
			data_mod=numpy.loadtxt('SNcc/chieffi/6E-3/'+str(el)+'.txt') #Open file
		elif (input_model == "Chieffi_2E-2"):
			data_mod=numpy.loadtxt('SNcc/chieffi/2E-2/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0.001"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0.004"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0.008"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0.008/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0.02"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNcc_0.05"):
			data_mod=numpy.loadtxt('SNcc/N13_SNcc/0.05/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_PISNe_0"):
			data_mod=numpy.loadtxt('SNcc/N13_PISNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "N13_SNe_0"):
			data_mod=numpy.loadtxt('SNcc/N13_SNe/0/'+str(el)+'.txt') #Open file
		else:
			print "Error: The model", input_model, "does not exist, or could not be found."
			sys.exit()

		yiel=numpy.arange(len(m))
		yiel=0.0*yiel
		if len(data_mod.shape)==1:
			yiel = data_mod
		else:
			for j in range(len(data_mod[:,0])):       #Merge the isotopes of a same element
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
def checkFe(cl_abun):
	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")
	if cl_abun[Fe_index,1] !=1.0:
		print "\n\nWarning: Fe abundance is not set to 1.0 in the imput file!"
		print "         Since the sum of the Fe from SNe models always scales to the Fe abundance in the imput file,"
		print "         The fit might be biased in case of Fe large uncertainties... "
		print "         It is strongly recommended to scale all the abundances to the Fe value!!\n\n"
	return 




## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).
def fitfunc1(modmatrix, cl_abun):
	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	SNe=1.0/modmatrix[Fe_index,0]
	return SNe



## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).
def fitfunc2(SNe,modmatrix, cl_abun):
	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	SNe[1]=(cl_abun[Fe_index,1]-SNe[0]*modmatrix[0,Fe_index])/modmatrix[1,Fe_index]
	return (cl_abun[:,1] - numpy.dot(SNe,modmatrix))




## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).

def fitfunc3(SNe,modmatrix, cl_abun):
	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	SNe[2]=(cl_abun[Fe_index,1]-SNe[0]*modmatrix[0,Fe_index]-SNe[1]*modmatrix[1,Fe_index])/modmatrix[2,Fe_index]
	return (cl_abun[:,1] - numpy.dot(SNe,modmatrix))





## Function used to fit the data to the models. This defines a
## set of linear equations (by multiplying the matrices of parameters
## (SNe) and of the predicted yields (modmatrix)).

def fitfunc4(SNe,modmatrix, cl_abun):
	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	SNe[3]=(cl_abun[Fe_index,1]-SNe[0]*modmatrix[0,Fe_index]-SNe[1]*modmatrix[1,Fe_index]-SNe[2]*modmatrix[2,Fe_index])/modmatrix[3,Fe_index]
	return (cl_abun[:,1] - numpy.dot(SNe,modmatrix))





## Removes the '\n' at the end of each item in a list.
def remove_end_char(s):
	return s[:-1]



## Give to each considered element its respective name.
def elem_name(elements):
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
	return elem_name












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
