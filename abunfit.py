#!/usr/bin/python

import os
import sys
import math
import numpy
import scipy.optimize as optimization
import matplotlib
import matplotlib.pyplot as plt

def main():

	####### User parameters #######

	abun_table = "solar_tables/lodders09.txt"	# Table of reference (solar) abundances

	alpha_SNcc_equiv_AGB = True	# If true, AGB model(s) follow the same IMF as for SNcc models
	alpha_AGB = -2.35		# Slope of the IMF for AGB model(s), if alpha_SNcc_equiv_AGB = False

	default_color_SNIa = "#4881ea"	# In plotting mode, the default color for the first SNIa model (if any)
	default_color_SNcc = "violet"	# In plotting mode, the default color for the first SNcc model (if any)
	default_color_AGB = "gold"	# In plotting mode, the default color for the first AGB model (if any)

	###############################




	if len(sys.argv) < 4:
		sys.exit("Error: Uncorrect arguments.\n\n \
		PROGRAM: abunfit.py\n\n\n \
		SYNOPSIS: python abunfit.py [input_file] [number_of_models] [alpha] [SNe_model_1] [SNe_model_2] ... (--disable-plot) \n\n\n \
		DESCRIPTION: Fits a set of intra-cluster abundances with a combination of \n \
		supernovae (SNe) yields models. Several models (either SNcc or SNIa) are available \n \
		and can be updated or created.\n\n\n \
		ARGUMENTS: \n\n \
		        [input_file] -- (Location of the) Text input file of the set of measured  \n \
		                        ICM abundances. This file must contain 3 columns:  \n \
		                        (1) the element number (Z), \n \
		                        (2) the measured abundance, and \n \
		                        (3) its associated (symmetrical) uncertainties.\n\n \
		        [number_of_models] -- Integer indicating the number of models to\n \
 			                       be fitted simultaneously (min: 1, max: 4).\n\n \
		        [alpha] -- Slope index of the assumed initial mass function (IMF). \n \
		                   Ex: for Salpeter IMF, alpha = -2.35.\n\n \
		        [model<1/2/...>] -- Name of the model. Currently available models: \n \
		                            ==============\n \
		                            ==== SNIa ====\n \
		                            ==============\n \
		                            == Iwamoto et al. (1999)\n \
		                             - Iw99_W7\n \
		                             - Iw99_W70\n \
		                             - Iw99_WDD1\n \
		                             - Iw99_WDD2\n \
		                             - Iw99_WDD3\n \
		                             - Iw99_CDD1\n \
		                             - Iw99_CDD2\n \
		                            == Badenes et al. (2006)\n \
		                               (Based on Tycho's best spectral fits...)\n \
		                             - Ba06_DDTa\n \
		                             - Ba06_DDTb\n \
		                             - Ba06_DDTc\n \
		                             - Ba06_DDTd\n \
		                             - Ba06_DDTe\n \
		                             - Ba06_DDTf\n \
		                            == Maeda et al. (2010)\n \
		                               (2-D explosion models...)\n \
		                             - Ma10_C-DEF\n \
		                             - Ma10_C-DDT\n \
		                             - Ma10_O-DDT\n \
		                            == Seitenzahl et al. (2013)\n \
		                               (3-D delayed-detonation explosion models...)\n \
		                             - Se13_N1\n \
		                             - Se13_N3\n \
		                             - Se13_N5\n \
		                             - Se13_N10\n \
		                             - Se13_N20\n \
		                             - Se13_N40\n \
		                             - Se13_N100H\n \
		                             - Se13_N100\n \
		                             - Se13_N100L\n \
		                             - Se13_N150\n \
		                             - Se13_N200\n \
		                             - Se13_N300C\n \
		                             - Se13_N1600\n \
		                             - Se13_N1600C\n \
		                             - Se13_N100_Z0.5\n \
		                             - Se13_N100_Z0.1\n \
		                             - Se13_N100_Z0.01\n \
		                            == Thielmann et al. (2010) \n \
		                            - Thielmann_03 \n \
		                            == Fink et al. (2014)\n \
		                               (3-D deflagration explosion models...)\n \
		                             - Fi14_N1def\n \
		                             - Fi14_N3def\n \
		                             - Fi14_N5def\n \
		                             - Fi14_N10def\n \
		                             - Fi14_N20def\n \
		                             - Fi14_N40def\n \
		                             - Fi14_N100Hdef\n \
		                             - Fi14_N100def\n \
		                             - Fi14_N100Ldef\n \
		                             - Fi14_N150def\n \
		                             - Fi14_N200def\n \
		                             - Fi14_N300Cdef\n \
		                             - Fi14_N1600def\n \
		                             - Fi14_N1600Cdef\n \
		                            == Ohlmann et al. (2014)\n \
		                               (3-D delayed-detonation expl. models, varying C fraction...)\n \
		                             - Oh14_DDT8_N100_c50\n \
		                             - Oh14_DDT8_N100_rpc20\n \
		                             - Oh14_DDT8_N100_rpc32\n \
		                             - Oh14_DDT8_N100_rpc40\n \
		                            == Leung & Nomoto (2018)\n \
		                               (2-D explosion models with varying initial metallicities...)\n \
		                             - Le18_050-1-c3-1P\n \
		                             - Le18_100-1-c3-1P\n \
		                             - Le18_100-0-c3\n \
		                             - Le18_100-0.1-c3\n \
		                             - Le18_100-0.5-c3\n \
		                             - Le18_100-1-c3\n \
		                             - Le18_100-2-c3\n \
		                             - Le18_100-3-c3\n \
		                             - Le18_100-5-c3\n \
		                             - Le18_300-1-c3-1P\n \
		                             - Le18_300-0-c3\n \
		                             - Le18_300-0.1-c3\n \
		                             - Le18_300-0.5-c3\n \
		                             - Le18_300-1-c3\n \
		                             - Le18_300-2-c3\n \
		                             - Le18_300-3-c3\n \
		                             - Le18_300-5-c3\n \
		                             - Le18_500-1-c3-1P\n \
		                             - Le18_500-0-c3\n \
		                             - Le18_500-0.1-c3\n \
		                             - Le18_500-0.5-c3\n \
		                             - Le18_500-1-c3\n \
		                             - Le18_500-2-c3\n \
		                             - Le18_500-3-c3\n \
		                             - Le18_500-5-c3\n \
		                            ==== Double-detonation (He shell) ====\n \
		                            == Sim et al. (2012)\n \
		                             - Si12_CSDD-L\n \
		                             - Si12_CSDD-S\n \
		                             - Si12_ELDD-L\n \
		                             - Si12_ELDD-S\n \
		                             - Si12_HeD-L\n \
		                             - Si12_HeD-S\n \
		                            ==== Gravitationally confined detonation ====\n \
		                            == Seitenzahl et al. (2016)\n \
		                             - Se16_GCD200\n \
		                            ==== Oxygen-neon WD ====\n \
		                            == Marquardt et al. (2015)\n \
		                             - Ma15_CO15e7\n \
		                             - Ma15_ONe10e7\n \
		                             - Ma15_ONe13e7\n \
		                             - Ma15_ONe15e7\n \
		                             - Ma15_ONe17e7\n \
		                             - Ma15_ONe20e7\n \
		                            ==== Hybrid WD (CO and ONe layers) ====\n \
		                            == Kromer et al. (2015)\n \
		                             - Kromer_N5_hybrid\n \
		                            ==== Ca-rich gap transient SNe ====\n \
		                            == Waldman et al. (2011)\n \
		                             - Wa11_CO.45HE.2\n \
		                             - Wa11_CO.55HE.2\n \
		                             - Wa11_CO.5HE.15\n \
		                             - Wa11_CO.5HE.2\n \
		                             - Wa11_CO.5HE.2C.3\n \
		                             - Wa11_CO.5HE.2N.02\n \
		                             - Wa11_CO.5HE.3\n \
		                             - Wa11_CO.6HE.2\n \
		                            == Sim et al. (2010)\n \
		                               (Detonation, single deg. channel...)\n \
		                             - Si10_det_0.81\n \
		                             - Si10_det_0.88\n \
		                             - Si10_det_0.97\n \
		                             - Si10_det_1.06\n \
		                             - Si10_det_1.06_0.075Ne\n \
		                             - Si10_det_1.15\n \
		                            ==== Sub-Chandrasekhar DD merger ====\n \
		                            == Pakmor et al. (2010)\n \
		                               (Violent merger of 0.9+0.9 M_sun...)\n \
		                             - Pa10_09_09\n \
		                            == Pakmor et al. (2012)\n \
		                               (Violent merger of 1.1+0.9 M_sun...)\n \
		                             - Pa12_11_09\n \
		                            == Kromer et al. (2013)\n \
		                               (Violent merger of 0.9+0.76 M_sun...)\n \
		                             - Kr13_09_076\n \
		                            == Kromer et al. (2013)\n \
		                               (Same, with lower initial metallicity...)\n \
		                             - Kr13_09_076_Z0.01\n \
		                            == Shen et al. (2018)\n \
		                               (Dynamically-driven double-degenerate double-detonation)\n \
		                             - Sh18_M[08/085/09/10/11]_[3070/5050]_Z[0/0005/001/002]_[01/1]\n \
		                            ==============\n \
		                            ==== SNcc ====\n \
		                            ==============\n \
		                            == Chieffi & Limongi (2004)\n \
		                               (One initial metallicity per model...)\n \
		                             - Ch04_0\n \
		                             - Ch04_1E-6\n \
		                             - Ch04_1E-4\n \
		                             - Ch04_1E-3\n \
		                             - Ch04_6E-3\n \
		                             - Ch04_2E-2\n \
		                            == Nomoto et al. (2006)\n \
		                               (One initial metallicity per model...)\n \
		                             - No06_0\n \
		                             - No06_0.001\n \
		                             - No06_0.004\n \
		                             - No06_0.02\n \
		                            == Romano et al. (2010) \n \
		                             - Ro10_0 \n \
		                             - Ro10_0.000002 \n \
		                             - Ro10_0.0002 \n \
		                             - Ro10_0.002 \n \
		                             - Ro10_0.02 \n \
		                            == Nomoto et al. (2013)\n \
		                               (One initial metallicity per model...)\n \
		                             = SNcc: 11 - 40 (140) M_sun\n \
		                              - No13_SNcc_0\n \
		                              - No13_SNcc_0.001\n \
		                              - No13_SNcc_0.004\n \
		                              - No13_SNcc_0.008\n \
		                              - No13_SNcc_0.02\n \
		                              - No13_SNcc_0.05\n \
		                             = PISNe: 140 - 300 M_sun\n \
		                              - No13_PISNe_0\n \
		                             = All SNe (SNcc+PISNe): 11 - 300 M_sun\n \
		                              - No13_SNe_0\n \
		                             = HNe: 20 - 40 (140) M_sun\n \
		                              - No13_HNe_0\n \
		                              - No13_HNe_0.001\n \
		                              - No13_HNe_0.004\n \
		                              - No13_HNe_0.008\n \
		                              - No13_HNe_0.02\n \
		                              - No13_HNe_0.05\n \
		                            == Heger & Woosley (2002,2010)\n \
		                               (Initial metallicity always 0!)\n \
		                             = SNcc: 10 - 100 M_sun\n \
		                              - He0210_SNcc_0\n \
		                             = PISNe: 140 - 260 M_sun\n \
		                              - He0210_PISNe_0\n \
		                             = All SNe (SNcc+PISNe): 10 - 260 M_sun\n \
		                              - He0210_SNe_0\n \
		                            == Sukhbold et al. (2016)\n \
		                              - Su16_N20\n \
		                              - Su16_W18\n \
		                            ===================\n \
		                            ==== AGB stars ====\n \
		                            ===================\n \
		                            == Karakas et al. (2010) \n \
		                             - K10_AGB_0.0001 \n \
		                             - K10_AGB_0.004 \n \
		                             - K10_AGB_0.008 \n \
		                             - K10_AGB_0.02 \n \
		                            == Nomoto et al. (2013)\n \
		                               (One initial metallicity per model...)\n \
		                             = AGB: 0.9 - 6.5 M_sun (adapted from Karakas 2010)\n \
		                              - No13_AGB_0.001\n \
		                              - No13_AGB_0.004\n \
		                              - No13_AGB_0.008\n \
		                              - No13_AGB_0.02\n\n \
		        --disable-plot -- If this argument ends the command line, \n \
		        		  no plot is displayed after the fit \n \
		        		  (useful for fitting in series, e.g. using serial_fitter.sh) \n\n\n \
		EXAMPLES: python abunfit.py data/CHEERS_Mernier18b.dat 2 -2.35 No13_SNcc_0.02 Se13_N100\n \
			  	--> Fits two models (1 SNcc + 1 SNIa) to the CHEERS abundance ratios \n \
			  	    (Mernier et al. 2018) assuming a Salpeter IMF. \n\n \
			  python abunfit.py data/CHEERS_Mernier18b.dat 2 -2.35 No13_SNcc_0.02 Se13_N100 --disable-plot\n \
			  	--> Same as above with the plotting option disabled. \n\n \
			  python abunfit.py data/proto-solar_Lodders09.dat 3 -2.35 No13_AGB_0.02 No13_SNcc_0.02 Se13_N100\n \
			  	--> Fits three models (1 AGB + 1 SNcc + 1 SNIa) to the proto-solar ratios \n \
			  	    (Lodders et al. 2009) assuming a Salpeter IMF. \n\n \
			  python abunfit.py data/proto-solar_Lodders09.dat 3 -2.35 No13_AGB_0.02 No13_SNcc_0.02 Se13_N100\n \
			  	--> Same as above, now assuming a top-heavy IMF. \n\n\n \
		VERSION: 3.0 (October 2023)\n \
		         2.0 (November 2018)\n \
		         1.5 (June 2018)\n \
		         1.4 (January 2018)\n \
		         1.3 (December 2017)\n \
		         1.2 (February 2016)\n \
		         1.1 (January 2016)\n \
		         1.0 (July 2015)\n\n\n \
		AUTHOR: Francois Mernier\n\n\n")


	plot_arg = "t"

	inputfile=sys.argv[1]
	number_of_models=sys.argv[2]
	number_of_models = int(number_of_models)
	alpha=sys.argv[3]
	input_mod_1=sys.argv[4]
	if number_of_models==1:
		if len(sys.argv) > 5:
			plot_arg=sys.argv[5]
	if number_of_models>1:
		input_mod_2=sys.argv[5]
		if len(sys.argv) > 6:
			plot_arg=sys.argv[6]
	if number_of_models>2:
		input_mod_3=sys.argv[6]
		if len(sys.argv) > 7:
			plot_arg=sys.argv[7]
	if number_of_models>3:
		input_mod_4=sys.argv[7]
		if len(sys.argv) > 8:
			plot_arg=sys.argv[8]
	if number_of_models>4 or number_of_models<1:
		sys.exit("Error: Abunfit can fit only 1, 2, 3 or 4 models simultaneously.")




	# Setting the plot mode ON or OFF (based on whether the plot_arg arggument exists and is --disable-plot
	if plot_arg=="--disable-plot":
		plotting_mode = False # Display numbers on the terminal only. Useful to fit in series (e.g. serial_fitter.sh)
	else:
		plotting_mode = True # Plot a figure with the best-fit models and save it into "fig_output.pdf"





##############################################################################################
#	Load files and tables...
##############################################################################################



	# Read cluster measurements data.
	f1 = open(inputfile, 'r')
	cl_abun = numpy.loadtxt(f1)

	elements = cl_abun[:,0].astype(int)
	data = cl_abun[:,1]
	errors = cl_abun[:,2]

	# IMF assumption: Salpeter
	#alpha = -2.35
	alpha = float(alpha)
	if alpha_SNcc_equiv_AGB == True:
		alpha_AGB = alpha



	# Give to each considered element its respective name.
	cl_elem_name = elem_name(elements)




	# Read atomic masses. 
	# mass[i] = atomic mass of the ith element.
	f2 = open("atomic_masses.txt", 'r')
	atom_mass = numpy.loadtxt(f2)
	mass = atom_mass.T[1][elements-1]






	# Read abundance tables (solarref et al. 2009). 
	# solarref[i] = protosolar abun. of the ith element (rel. to Si==10^6).
	f3 = open(abun_table, 'r')
	abun_file = numpy.loadtxt(f3)
	solarref = abun_file.T[1][elements-1]











##############################################################################################
#	Set up (and, if relevant, integrate) yield models...
##############################################################################################




	# Read list of available SNe models. 
	list_models_SNIa = numpy.loadtxt("list_models_SNIa.txt", dtype="str").tolist()
	list_models_SNcc = numpy.loadtxt("list_models_SNcc.txt", dtype="str").tolist()
	list_models_AGB = numpy.loadtxt("list_models_AGB.txt", dtype="str").tolist()



	# Load SNIa and/or SNcc models...

	if (input_mod_1 in list_models_SNIa):
		model1 = sum_SNIa(input_mod_1,elements)
	elif (input_mod_1 in list_models_SNcc):
		model1 = integrate_SNcc(input_mod_1,alpha,cl_elem_name)
	elif (input_mod_1 in list_models_AGB):
		model1 = integrate_AGB(input_mod_1,alpha_AGB,cl_elem_name)
	else:
		print("Error: The model", input_mod_1, "does not exist, or could not be found.")
		sys.exit()


	if number_of_models>1:
		if (input_mod_2 in list_models_SNIa):
			model2 = sum_SNIa(input_mod_2,elements)
		elif (input_mod_2 in list_models_SNcc):
			model2 = integrate_SNcc(input_mod_2,alpha,cl_elem_name)
		elif (input_mod_2 in list_models_AGB):
			model2 = integrate_AGB(input_mod_2,alpha_AGB,cl_elem_name)
		else:
			print("Error: The model", input_mod_2, "does not exist, or could not be found.")
			sys.exit()


	if number_of_models>2:
		if (input_mod_3 in list_models_SNIa):
			model3 = sum_SNIa(input_mod_3,elements)
		elif (input_mod_3 in list_models_SNcc):
			model3 = integrate_SNcc(input_mod_3,alpha,cl_elem_name)
		elif (input_mod_3 in list_models_AGB):
			model3 = integrate_AGB(input_mod_3,alpha_AGB,cl_elem_name)
		else:
			print("Error: The model", input_mod_3, "does not exist, or could not be found.")
			sys.exit()


	if number_of_models>3:
		if (input_mod_4 in list_models_SNIa):
			model4 = sum_SNIa(input_mod_4,elements)
		elif (input_mod_4 in list_models_SNcc):
			model4 = integrate_SNcc(input_mod_4,alpha,cl_elem_name)
		elif (input_mod_4 in list_models_AGB):
			model4 = integrate_AGB(input_mod_4,alpha_AGB,cl_elem_name)
		else:
			print("Error: The model", input_mod_4, "does not exist, or could not be found.")
			sys.exit()

	

	print("Model: ", input_mod_1)
	print(model1)
	if number_of_models>1:
		print("Model: ", input_mod_2)
		print(model2)
	if number_of_models>2:
		print("Model: ", input_mod_3)
		print(model3)
	if number_of_models>3:
		print("Model: ", input_mod_4)
		print(model4)



	x=numpy.zeros(len(elements))
	if number_of_models>1:
		y=numpy.zeros(len(elements))
	if number_of_models>2:
		z=numpy.zeros(len(elements))
	if number_of_models>3:
		w=numpy.zeros(len(elements))


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=model1[i]/(mass[i]*solarref[i])
		if number_of_models>1:
			y[i]=model2[i]/(mass[i]*solarref[i])
		if number_of_models>2:
			z[i]=model3[i]/(mass[i]*solarref[i])
		if number_of_models>3:
			w[i]=model4[i]/(mass[i]*solarref[i])


	









##############################################################################################
#	Fit...
##############################################################################################



	SNe=numpy.arange(number_of_models)
	SNe=SNe+1000
	x0=SNe

	if number_of_models==1:
		modmatrix = numpy.vstack((x))
		#print optimization.leastsq(fitfunc1, x0, args=(modmatrix, data))
		bestSNe=numpy.array([0])
		bestSNe[0]=fitfunc1(modmatrix, cl_abun)
		bestfit=numpy.dot(modmatrix,bestSNe)

	if number_of_models==2:
		modmatrix = numpy.vstack((x,y))
		#print optimization.leastsq(fitfunc2, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc2, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)

	if number_of_models==3:
		modmatrix = numpy.vstack((x,y,z))
		#print optimization.leastsq(fitfunc3, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc3, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)

	if number_of_models==4:
		modmatrix = numpy.vstack((x,y,z,w))
		#print optimization.leastsq(fitfunc4, x0, args=(modmatrix, data))
		bestSNe=optimization.leastsq(fitfunc4, x0, args=(modmatrix, cl_abun))[0]
		bestfit=numpy.dot(bestSNe,modmatrix)


	if number_of_models==2:
		sxz=0.
		syz=0.
		sxx=0.
		syy=0.
		sxy=0.

		for i in range(len(elements)):
			#sxz = sxz + x[i]*data[i]/(errors[i]**2)
			#syz = syz + y[i]*data[i]/(errors[i]**2)
			sxx = sxx + x[i]**2/(errors[i]**2)
			syy = syy + y[i]**2/(errors[i]**2)
			sxy = sxy + x[i]*y[i]/(errors[i]**2)
	
		#a = (sxz*syy - syz*sxy)/(syy*sxx - sxy**2)
		#b = (sxx*syz - sxz*sxy)/(syy*sxx - sxy**2)
    
		delta = syy*sxx - sxy**2
		siga2 = syy / delta
		sigb2 = sxx / delta

		#ratio = a / (a + b)
		errn = math.sqrt(siga2 + sigb2)
		errr = math.sqrt((errn/(bestSNe[0]+bestSNe[1]))**2 + (siga2/(bestSNe[0]**2)))








	#Check for possible (unphysical) negative contribution...
	for i in range(len(bestSNe)):
		if bestSNe[i] < 0.0:
			print("Error: One model has a negative contribution to the enrichment. This is physically not allowed.")
			chi = 9999999
			chiout = open("currentchi.txt", 'w+')
			chiout.write("%s" % int(chi*100000))
			sys.exit()




	#Compute chi2...
	chi = 0.0
	print("Optimal parameters:", bestSNe)
	if number_of_models==2:
		print("Ratio:", bestSNe[0] / (bestSNe[0]+bestSNe[1]), "+/-", errr)
	for i in range(len(elements)):
		if number_of_models==1:
			chi += (data[i] - bestSNe[0]*x[i])**2 / errors[i]**2
		if number_of_models==2:
			chi += (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i])**2 / errors[i]**2
		if number_of_models==3:
			chi += (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i])**2 / errors[i]**2
		if number_of_models==4:
			chi += (data[i] - bestSNe[0]*x[i] - bestSNe[1]*y[i] - bestSNe[2]*z[i] - bestSNe[3]*w[i])**2 / errors[i]**2

	print("chi^2 / d.o.f. = ", float(chi),'/',(len(elements)- number_of_models), "=", \
	       chi / (len(elements) - number_of_models))

	chiout = open("currentchi.txt", 'w+') 
	chiout.write("%s" % int(chi*100000)) # Used for serial_fitter.sh. *100000-> to distinguish the decimals (will be corrected in sort_fitting_results.py)

	chiout.close()



##############################################################################################
#	Print results in a QDP/text file (output.qdp)...
##############################################################################################


	print("--------------------------QDP FILE--------------------------")
	print("skip single")
	print("line on 2,3,4,5,6")
	print("mark 0 on 1")
	print("r y -0.1 3")
	print("la x Atomic Number")
	print("la y Abundance (proto-solar)")
	print("READ Serr 1,2")
	print("! Abundance measurements")
	for el in range(len(elements)):
		print(elements[el], 0, data[el], errors[el])
	print("NO")
	print("! Total of models")
	for el in range(len(elements)):
		print(elements[el], 0, bestfit[el], 0)
	print("NO")
	print("! Model:", input_mod_1)
	for el in range(len(elements)):
		print(elements[el], 0, bestSNe[0]*x[el], 0)
	if number_of_models>1:
		print("NO")
		print("! Model:", input_mod_2)
		for el in range(len(elements)):
			print(elements[el], 0, bestSNe[1]*y[el], 0)
	if number_of_models>2:
		print("NO")
		print("! Model:", input_mod_3)
		for el in range(len(elements)):
			print(elements[el], 0, bestSNe[2]*z[el], 0)
	if number_of_models>3:
		print("NO")
		print("! Model:", input_mod_4)
		for el in range(len(elements)):
			print(elements[el], 0, bestSNe[3]*w[el], 0)





	#Write into QDP file...
	qdpout = open("output.qdp", 'w+')
	qdpout.write("skip single\n")
	qdpout.write("line on 2,3,4,5,6,7\n")
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









##############################################################################################
#	Plotting with matplotlib...
##############################################################################################


	if plotting_mode == True:
		fig = plt.figure(figsize=(8, 5))
		fig.subplots_adjust(bottom=0.12, right=0.95, top=0.95)

		xmax = len(elements)+1
		ymax = rescale_ymax(number_of_models, bestfit, data, errors)

		width_bars=0.60
		plt.axis([0.5, xmax, 0.0, ymax])
		#plt.xlabel('Atomic numer', fontsize=14)
		plt.ylabel("X/Fe Abundance ratio (proto-solar)", fontsize=12)
		plt.tick_params(labelsize=13)
		xaxis = numpy.arange(len(elements)+1)
		xaxis = numpy.delete(xaxis, 0)
		xlabels=cl_elem_name
		plt.xticks(xaxis+width_bars/2., xlabels)

		fraction_label = label_SNfrac(input_mod_1, input_mod_2, list_models_SNIa, list_models_SNcc, list_models_AGB)


		# Plot histograms and data points (depending on the number of models)

		if number_of_models==1:
			color1 = set_colors_1(input_mod_1, default_color_SNIa, default_color_SNcc, default_color_AGB)
			p1 = plt.bar(tuple(xaxis), tuple(bestSNe[0]*x[:]), width_bars, align="edge", color=color1, zorder=1)
			ptot = plt.bar(tuple(xaxis), tuple(bestfit[:]), width_bars, align="edge", color='none', edgecolor='black', zorder=2)
			input_mod_1 = rename_model(input_mod_1, list_models_SNIa, list_models_SNcc, list_models_AGB)


		if number_of_models==2:
			color2 = set_colors_2(input_mod_1, input_mod_2, default_color_SNIa, default_color_SNcc, default_color_AGB)
			p1 = plt.bar(tuple(xaxis), tuple(bestSNe[0]*x[:]), width_bars, align="edge", color=color2[0], zorder=1)
			p2 = plt.bar(tuple(xaxis), tuple(bestSNe[1]*y[:]), width_bars, align="edge", color=color2[1], bottom=bestSNe[0]*x[:], zorder=1)
			input_mod_1 = rename_model(input_mod_1, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_2 = rename_model(input_mod_2, list_models_SNIa, list_models_SNcc, list_models_AGB)
			ptot = plt.bar(tuple(xaxis), tuple(bestfit[:]), width_bars, align="edge", color='none', edgecolor='black', zorder=2)

		if number_of_models==3:
			color3 = set_colors_3(input_mod_1, input_mod_2, input_mod_3, default_color_SNIa, default_color_SNcc, default_color_AGB)
			p1 = plt.bar(tuple(xaxis), tuple(bestSNe[0]*x[:]), width_bars, align="edge", color=color3[0], zorder=1)
			p2 = plt.bar(tuple(xaxis), tuple(bestSNe[1]*y[:]), width_bars, align="edge", color=color3[1], bottom=bestSNe[0]*x[:], zorder=1)
			p3 = plt.bar(tuple(xaxis), tuple(bestSNe[2]*z[:]), width_bars, align="edge", color=color3[2], bottom=bestSNe[0]*x[:]+bestSNe[1]*y[:], zorder=1)
			input_mod_1 = rename_model(input_mod_1, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_2 = rename_model(input_mod_2, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_3 = rename_model(input_mod_3, list_models_SNIa, list_models_SNcc, list_models_AGB)
			ptot = plt.bar(tuple(xaxis), tuple(bestfit[:]), width_bars, align="edge", color='none', edgecolor='black', zorder=2)

		if number_of_models==4:
			color4 = set_colors_4(input_mod_1, input_mod_2, input_mod_3, input_mod_4, default_color_SNIa, default_color_SNcc, default_color_AGB)
			print(bestSNe[3]*w[:])
			print(color4)
			p1 = plt.bar(tuple(xaxis), tuple(bestSNe[0]*x[:]), width_bars, align="edge", color=color4[0], zorder=1)
			p2 = plt.bar(tuple(xaxis), tuple(bestSNe[1]*y[:]), width_bars, align="edge", color=color4[1], bottom=bestSNe[0]*x[:], zorder=1)
			p3 = plt.bar(tuple(xaxis), tuple(bestSNe[2]*z[:]), width_bars, align="edge", color=color4[2], bottom=bestSNe[0]*x[:]+bestSNe[1]*y[:], zorder=1)
			p4 = plt.bar(tuple(xaxis), tuple(bestSNe[3]*w[:]), width_bars, align="edge", color=color4[3], bottom=bestSNe[0]*x[:]+bestSNe[1]*y[:]+bestSNe[2]*z[:], zorder=1)
			input_mod_1 = rename_model(input_mod_1, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_2 = rename_model(input_mod_2, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_3 = rename_model(input_mod_3, list_models_SNIa, list_models_SNcc, list_models_AGB)
			input_mod_4 = rename_model(input_mod_4, list_models_SNIa, list_models_SNcc, list_models_AGB)
			ptot = plt.bar(tuple(xaxis), tuple(bestfit[:]), width_bars, align="edge", color='none', edgecolor='black', zorder=2)



		item1 = plt.errorbar(xaxis[:]+width_bars/2., data[:], xerr=0.0, yerr=errors[:], ls=' ', linewidth=1.0, color='black', markersize=6, fmt='s', zorder=3)
		plt.axhline(y=1.0, xmin=0, xmax=30, color='black', linestyle='--', linewidth=1.0, alpha=0.4, zorder=0)
		plt.annotate(r'$\chi^2$/d.o.f. = '+str(numpy.round_(float(chi),1))+"/"+str((len(elements)- number_of_models)), xy=(0.95*xmax, 0.9*ymax), horizontalalignment='right', color='black', fontsize=14)

		plt.annotate(fraction_label+str(round(bestSNe[0] / (bestSNe[0]+bestSNe[1]),4)), xy=(0.9*xmax, 0.8*ymax), horizontalalignment='right', color='black', fontsize=16)




		# Plot legend (depending on the number of models)

		if number_of_models==1:
			plt.legend([item1, p1], ["Data", input_mod_1], ncol=1, numpoints=1, loc=(0.02,0.85), fontsize=10)
		if number_of_models==2:
			plt.legend([item1, p1, p2], ["Data", input_mod_1, input_mod_2], ncol=1, numpoints=1, loc=(0.02,0.80), fontsize=10)
		if number_of_models==3:
			plt.legend([item1, p1, p2, p3], ["Data", input_mod_1, input_mod_2, input_mod_3], ncol=1, numpoints=1, loc=(0.02,0.75), fontsize=10)
		if number_of_models==4:
			plt.legend([item1, p1, p2, p3, p4], ["Data", input_mod_1, input_mod_2, input_mod_3, input_mod_4], ncol=1, numpoints=1, loc=(0.02,0.70), fontsize=10)

		plt.show()
		fig.savefig('fig_output.pdf', format='pdf')



























	#Last check, to be sure that Fe is normalized to 1.0 ...
	checkFe(cl_abun)























##############################################################################################
#	ADDITIONAL FUNCTIONS
##############################################################################################



## Integrates the SNcc products over a mass range of 13-40 M_sun
## and assuming a certain IMF.
def integrate_SNcc( input_model, alpha , elem_name ):

	if (input_model == "No06_0") or (input_model == "No06_0.001") or \
		(input_model == "No06_0.004") or (input_model == "No06_0.02") or \
		(input_model == "No06_0_undecayed") or (input_model == "No06_0.001_undecayed") or \
		(input_model == "No06_0.004_undecayed") or (input_model == "No06_0.02_undecayed"):
		m=numpy.array([13,15,18,20,25,30,40])
	elif (input_model == "Ch04_0") or (input_model == "Ch04_1E-6") or \
		(input_model == "Ch04_1E-4") or (input_model == "Ch04_1E-3") or \
		(input_model == "Ch04_6E-3") or (input_model == "Ch04_2E-2"):
		m=numpy.array([13,15,20,25,30,35])
	elif (input_model == "No13_SNcc_0"):
		m=numpy.array([11,13,15,18,20,25,30,40,100,140])
	elif (input_model == "No13_SNcc_0.001") or (input_model == "No13_SNcc_0.004") or \
		(input_model == "No13_SNcc_0.008") or (input_model == "No13_SNcc_0.02") or \
		(input_model == "No13_SNcc_0.05"):
		m=numpy.array([13,15,18,20,25,30,40])
	elif (input_model == "No13_PISNe_0"):
		m=numpy.array([140,150,170,200,270,300])
	elif (input_model == "No13_SNe_0"):
		m=numpy.array([11,13,15,18,20,25,30,40,100,140,150,170,200,270,300])
	elif (input_model == "No13_HNe_0"):
		m=numpy.array([20,25,30,40,100,140])
	elif (input_model == "No13_HNe_0.001") or (input_model == "No13_HNe_0.004") or \
		(input_model == "No13_HNe_0.008") or (input_model == "No13_HNe_0.02") or \
		(input_model == "No13_HNe_0.05"):
		m=numpy.array([20,25,30,40])
	elif (input_model == "He0210_SNcc_0"):
		m=numpy.array([10,12,15,20,25,35,50,75,100])
	elif (input_model == "He0210_PISNe_0"):
		m=numpy.array([140,150,158,168,177,186,195,205,214,223,232,242,251,260])
	elif (input_model == "He0210_SNe_0"):
		m=numpy.array([10,12,15,20,25,35,50,75,100,140,150,158,168,177,186,195,205,214,223,232,242,251,260])
	elif (input_model == "Su16_N20"):
		m=numpy.array([12.25, 12.5, 12.75, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 
		14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.2, 15.7, 15.8, 
		15.9, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 
		17.3, 17.4, 17.5, 17.6, 17.7, 17.9, 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.7, 
		18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.7, 19.8, 20.1, 20.2, 20.3, 20.4, 
		20.5, 20.6, 20.8, 21.0, 21.1, 21.2, 21.5, 21.6, 21.7, 25.2, 25.3, 25.4, 25.5, 
		25.6, 25.7, 25.8, 25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 
		26.9, 27.0, 27.1, 27.2, 27.3, 27.4, 29.0, 29.1, 29.2, 29.6, 60, 80, 100, 120])
	elif (input_model == "Su16_W18"):
		m=numpy.array([12.25, 12.5, 12.75, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 
		14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.2, 15.7, 15.8, 
		16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.3, 
		17.4, 17.5, 17.6, 17.9, 18.1, 18.2, 18.3, 18.4, 18.5, 19.2, 19.3, 19.7, 19.8, 
		20.1, 20.2, 20.3, 20.4, 20.5, 20.8, 21.0, 21.1, 21.2, 21.5, 21.6, 25.2, 25.4, 
		25.5, 25.6, 25.7, 25.8, 25.9, 26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 27.0, 27.1, 
		27.2, 27.3, 60, 120])
	elif (input_model == "Ro10_0") or (input_model == "Ro10_0.000002") or (input_model == "Ro10_0.0002") or \
		(input_model == "Ro10_0.002") or (input_model == "Ro10_0.02"):
		m=numpy.array([1.])

	y=numpy.zeros(len(elem_name))

	alpha = float(alpha)
	#alpha=-2.35 (Salpeter)
	#alpha=-1.0 (Top)
    
  #for z in metal:     
	n=0

	for el in elem_name:
		if (input_model == "No06_0"):
			data_mod=numpy.loadtxt('SNcc/No06/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No06_0.001"):
			data_mod=numpy.loadtxt('SNcc/No06/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "No06_0.004"):
			data_mod=numpy.loadtxt('SNcc/No06/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "No06_0.02"):
			data_mod=numpy.loadtxt('SNcc/No06/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_0"):
			data_mod=numpy.loadtxt('SNcc/Ch04/0/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_1E-6"):
			data_mod=numpy.loadtxt('SNcc/Ch04/1E-6/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_1E-4"):
			data_mod=numpy.loadtxt('SNcc/Ch04/1E-4/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_1E-3"):
			data_mod=numpy.loadtxt('SNcc/Ch04/1E-3/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_6E-3"):
			data_mod=numpy.loadtxt('SNcc/Ch04/6E-3/'+str(el)+'.txt') #Open file
		elif (input_model == "Ch04_2E-2"):
			data_mod=numpy.loadtxt('SNcc/Ch04/2E-2/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.001"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.004"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.008"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0.008/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.02"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.05"):
			data_mod=numpy.loadtxt('SNcc/No13_SNcc/0.05/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_PISNe_0"):
			data_mod=numpy.loadtxt('SNcc/No13_PISNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNe_0"):
			data_mod=numpy.loadtxt('SNcc/No13_SNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_HNe_0"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_HNe_0.001"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_HNe_0.004"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_HNe_0.008"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0.008/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_HNe_0.02"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_SNcc_0.05"):
			data_mod=numpy.loadtxt('SNcc/No13_HNe/0.05/'+str(el)+'.txt') #Open file
		elif (input_model == "He0210_SNcc_0"):
			data_mod=numpy.loadtxt('SNcc/He0210_SNcc/0/'+str(el)+'.txt') #Open file
		elif (input_model == "He0210_PISNe_0"):
			data_mod=numpy.loadtxt('SNcc/He0210_PISNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "He0210_SNe_0"):
			data_mod=numpy.loadtxt('SNcc/He0210_SNe/0/'+str(el)+'.txt') #Open file
		elif (input_model == "Su16_N20"):
			data_mod=numpy.loadtxt('SNcc/Su16/N20/'+str(el)+'.txt') #Open file
		elif (input_model == "Su16_W18"):
			data_mod=numpy.loadtxt('SNcc/Su16/W18/'+str(el)+'.txt') #Open file
		elif (input_model == "Ro10_0"):
			data_mod=numpy.array([numpy.loadtxt('SNcc/Ro10/0/'+str(el)+'.txt')]) #Open file
		elif (input_model == "Ro10_0.000002"):
			data_mod=numpy.array([numpy.loadtxt('SNcc/Ro10/0.000002/'+str(el)+'.txt')]) #Open file
		elif (input_model == "Ro10_0.0002"):
			data_mod=numpy.array([numpy.loadtxt('SNcc/Ro10/0.0002/'+str(el)+'.txt')]) #Open file
		elif (input_model == "Ro10_0.002"):
			data_mod=numpy.array([numpy.loadtxt('SNcc/Ro10/0.002/'+str(el)+'.txt')]) #Open file
		elif (input_model == "Ro10_0.02"):
			data_mod=numpy.array([numpy.loadtxt('SNcc/Ro10/0.02/'+str(el)+'.txt')]) #Open file			
		else:
			print("Error: The model", input_model, "does not exist, or could not be found.")
			sys.exit()

		yiel=numpy.zeros(len(m))
		if len(data_mod.shape)==1:
			yiel=numpy.array(data_mod)
		else:
			yiel=numpy.array(numpy.sum(data_mod.T, axis=1))   

		sumt=0.
		sumn=0.

		if ("Ro10" in input_model):
			sumt=numpy.sum(yiel*m**alpha)
			sumn=numpy.sum(m**alpha)
		else:
			sumt=numpy.trapz(yiel*numpy.power(m,alpha), x=m)    #Integration "trapezoidal" (preferred)
			sumn=numpy.trapz(numpy.power(m,alpha), x=m)

		y[n]=sumt/sumn
		n=n+1


	return y









## Integrates the SNcc products over a mass range of 13-40 M_sun
## and assuming a certain IMF.
def integrate_AGB( input_model, alpha , elem_name ):

	if "K10_AGB" in input_model:
		m=numpy.array([1])
	else:
		m=numpy.array([0.9,1.0,1.25,1.5,1.75,1.9,2.0,2.25,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5])
	y=numpy.zeros(len(elem_name))

	alpha = float(alpha)
	#alpha=-2.35 (Salpeter)
	#alpha=-1.0 (Top)
    
  #for z in metal:     
	n=0

	for el in elem_name:
		if (input_model == "No13_AGB_0"):
			data_mod=numpy.loadtxt('AGB/No13_AGB/0/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_AGB_0.001"):
			data_mod=numpy.loadtxt('AGB/No13_AGB/0.001/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_AGB_0.004"):
			data_mod=numpy.loadtxt('AGB/No13_AGB/0.004/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_AGB_0.008"):
			data_mod=numpy.loadtxt('AGB/No13_AGB/0.008/'+str(el)+'.txt') #Open file
		elif (input_model == "No13_AGB_0.02"):
			data_mod=numpy.loadtxt('AGB/No13_AGB/0.02/'+str(el)+'.txt') #Open file
		elif (input_model == "K10_AGB_0.0001"):
			 data_mod=numpy.array([numpy.loadtxt('AGB/K10_AGB/0.0001/'+str(el)+'.txt')]) #Open file  
		elif (input_model == "K10_AGB_0.004"):
			 data_mod=numpy.array([numpy.loadtxt('AGB/K10_AGB/0.004/'+str(el)+'.txt')]) #Open file  
		elif (input_model == "K10_AGB_0.008"):
			 data_mod=numpy.array([numpy.loadtxt('AGB/K10_AGB/0.008/'+str(el)+'.txt')]) #Open file
		elif (input_model == "K10_AGB_0.02"):
			 data_mod=numpy.array([numpy.loadtxt('AGB/K10_AGB/0.02/'+str(el)+'.txt')]) #Open file	
		else:
			print("Error: The model", input_model, "does not exist, or could not be found.")
			sys.exit()

		yiel=numpy.zeros(len(m))
		if len(data_mod.shape)==1:
			yiel=numpy.array(data_mod)
		else:
			yiel=numpy.array(numpy.sum(data_mod.T, axis=1))   

		sumt=0.
		sumn=0.
	   	
		if ("K10_AGB" in input_model):
			sumt=numpy.sum(yiel*m**alpha)
			sumn=numpy.sum(m**alpha)
		else:
			sumt=numpy.trapz(yiel*numpy.power(m,alpha), x=m)    #Integration "trapezoidal" (preferred)
			sumn=numpy.trapz(numpy.power(m,alpha), x=m)

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














## Calculates the purely SNcc contribution to a fit from a AGB+SNcc component
## (letting the Fe relative contributions of the SNcc to 1)
def extract_SNcc_contrib_1model( input_AGBSNcc_model , alpha , elem_name , mass , solarref , elements , bestSNe , cl_abun ):

	if (input_AGBSNcc_model == "No13_AGB+SNcc_0"):
		input_SNcc_model = "No13_SNcc_0"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.001"):
		input_SNcc_model = "No13_SNcc_0.001"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.004"):
		input_SNcc_model = "No13_SNcc_0.004"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.008"):
		input_SNcc_model = "No13_SNcc_0.008"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.02"):
		input_SNcc_model = "No13_SNcc_0.02"
	else:
		print("Error: The model", input_AGBSNcc_model, "does not exist, or could not be found.")
		sys.exit()


	input_model1 = integrate_SNcc(input_SNcc_model,alpha,elem_name)



	x=numpy.arange(len(elements))
	x=x*0.0


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=input_model1[i]/(mass[i]*solarref[i])

	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	# Fit...
	SNe = 0.0

	modmatrix = numpy.vstack((x))
	SNe=1.0/modmatrix[Fe_index,0]

	return SNe*x














## Calculates the purely SNcc contribution to a fit from a AGB+SNcc component
## (letting the Fe relative contributions of the 2 models unchanged)
def extract_SNcc_contrib_2models( input_AGBSNcc_model, input_model2, alpha , elem_name , mass , solarref , elements , bestSNe , cl_abun ):

	if (input_AGBSNcc_model == "No13_AGB+SNcc_0"):
		input_SNcc_model = "No13_SNcc_0"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.001"):
		input_SNcc_model = "No13_SNcc_0.001"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.004"):
		input_SNcc_model = "No13_SNcc_0.004"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.008"):
		input_SNcc_model = "No13_SNcc_0.008"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.02"):
		input_SNcc_model = "No13_SNcc_0.02"
	else:
		print("Error: The model", input_AGBSNcc_model, "does not exist, or could not be found.")
		sys.exit()


	input_model1 = integrate_SNcc(input_SNcc_model,alpha,elem_name)



	x=numpy.arange(len(elements))
	x=x*0.0
	y=numpy.arange(len(elements))
	y=y*0.0


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=input_model1[i]/(mass[i]*solarref[i])
		y[i]=input_model2[i]/(mass[i]*solarref[i])

	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	# Fit...
	SNe=numpy.arange(2)
	SNe=SNe+1000
	SNe[1] = bestSNe[1]

	modmatrix = numpy.vstack((x,y))
	SNe[0]=(cl_abun[Fe_index,1]-bestSNe[1]*modmatrix[1,Fe_index])/modmatrix[0,Fe_index]

	return SNe[0]*x














## Calculates the purely SNcc contribution to a fit from a AGB+SNcc component
## (letting the Fe relative contributions of the 3 models unchanged)
def extract_SNcc_contrib_3models( input_AGBSNcc_model, input_model2, input_model3 , alpha , elem_name , mass , solarref , elements , bestSNe , cl_abun ):

	if (input_AGBSNcc_model == "No13_AGB+SNcc_0"):
		input_SNcc_model = "No13_SNcc_0"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.001"):
		input_SNcc_model = "No13_SNcc_0.001"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.004"):
		input_SNcc_model = "No13_SNcc_0.004"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.008"):
		input_SNcc_model = "No13_SNcc_0.008"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.02"):
		input_SNcc_model = "No13_SNcc_0.02"
	else:
		print("Error: The model", input_AGBSNcc_model, "does not exist, or could not be found.")
		sys.exit()


	input_model1 = integrate_SNcc(input_SNcc_model,alpha,elem_name)



	x=numpy.arange(len(elements))
	x=x*0.0
	y=numpy.arange(len(elements))
	y=y*0.0
	z=numpy.arange(len(elements))
	z=z*0.0


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=input_model1[i]/(mass[i]*solarref[i])
		y[i]=input_model2[i]/(mass[i]*solarref[i])
		z[i]=input_model3[i]/(mass[i]*solarref[i])

	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	# Fit...
	SNe=numpy.arange(3)
	SNe=SNe+1000
	SNe[1] = bestSNe[1]
	SNe[2] = bestSNe[2]

	modmatrix = numpy.vstack((x,y,z))
	SNe[0]=(cl_abun[Fe_index,1]-bestSNe[1]*modmatrix[1,Fe_index]-bestSNe[2]*modmatrix[2,Fe_index])/modmatrix[0,Fe_index]

	return SNe[0]*x














## Calculates the purely SNcc contribution to a fit from a AGB+SNcc component
## (letting the Fe relative contributions of the 4 models unchanged)
def extract_SNcc_contrib_4models( input_AGBSNcc_model , input_model2 , input_model3 , input_model4 , alpha , elem_name , mass , solarref , elements , bestSNe , cl_abun ):

	if (input_AGBSNcc_model == "No13_AGB+SNcc_0"):
		input_SNcc_model = "No13_SNcc_0"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.001"):
		input_SNcc_model = "No13_SNcc_0.001"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.004"):
		input_SNcc_model = "No13_SNcc_0.004"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.008"):
		input_SNcc_model = "No13_SNcc_0.008"
	elif (input_AGBSNcc_model == "No13_AGB+SNcc_0.02"):
		input_SNcc_model = "No13_SNcc_0.02"
	else:
		print("Error: The model", input_AGBSNcc_model, "does not exist, or could not be found.")
		sys.exit()


	input_model1 = integrate_SNcc(input_SNcc_model,alpha,elem_name)



	x=numpy.arange(len(elements))
	x=x*0.0
	y=numpy.arange(len(elements))
	y=y*0.0
	z=numpy.arange(len(elements))
	z=z*0.0
	w=numpy.arange(len(elements))
	w=w*0.0


	# Define X, Y (, Z, W): abun / supernova
	for i in range(len(elements)):
		x[i]=input_model1[i]/(mass[i]*solarref[i])
		y[i]=input_model2[i]/(mass[i]*solarref[i])
		z[i]=input_model3[i]/(mass[i]*solarref[i])
		w[i]=input_model4[i]/(mass[i]*solarref[i])

	Fe_index=-1
	for i in range(len(cl_abun[:,0])):	# Calculate the array index for Fe
		if cl_abun[i,0] ==26:
			Fe_index=i
	if Fe_index<1:
		sys.exit("Error: Fe abundance not found in", inputfile, ".")

	# Fit...
	SNe=numpy.arange(4)
	SNe=SNe+1000
	SNe[1] = bestSNe[1]
	SNe[2] = bestSNe[2]
	SNe[3] = bestSNe[3]

	modmatrix = numpy.vstack((x,y,z,w))
	SNe[0]=(cl_abun[Fe_index,1]-bestSNe[1]*modmatrix[1,Fe_index]-bestSNe[2]*modmatrix[2,Fe_index]-bestSNe[3]*modmatrix[3,Fe_index])/modmatrix[0,Fe_index]

	return SNe[0]*x














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
		print("\n\nWarning: Fe abundance is not set to 1.0 in the imput file!")
		print("         Since the sum of the Fe from SNe models always scales to the Fe abundance in the imput file,")
		print("         The fit might be biased in case of Fe large uncertainties... ")
		print("         It is strongly recommended to scale all the abundances to the Fe value!!\n\n")
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
	return ((cl_abun[:,1] - numpy.dot(SNe,modmatrix))/cl_abun[:,2])



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
	return ((cl_abun[:,1] - numpy.dot(SNe,modmatrix))/cl_abun[:,2])





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
	return ((cl_abun[:,1] - numpy.dot(SNe,modmatrix))/cl_abun[:,2])



## Removes the '\n' at the end of each item in a list.
def remove_end_char(s):
	return s[:-1]


## Gives to each considered element its respective name.
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



## Rescales the maximum value of the y-axis displayed in the final plot.
def rescale_ymax(number_of_models, bestfit, data, errors):
	bestfit1 = numpy.array_split(bestfit,2)[0]
	bestfit2 = numpy.array_split(bestfit,2)[1]
	data1 = numpy.array_split(data,2)[0]
	data2 = numpy.array_split(data,2)[1]
	errors1 = numpy.array_split(errors,2)[0]
	errors2 = numpy.array_split(errors,2)[1]
	if number_of_models==1:
		if numpy.max(bestfit1[:]) < 1.7 and numpy.max(data1[:]+errors1[:]) < 1.7 and \
		numpy.max(bestfit2[:]) < 1.75 and numpy.max(data2[:]+errors2[:]) < 1.75:
			ymax = 2.0
		elif numpy.max(bestfit1[:]) < 1.87 and numpy.max(data1[:]+errors1[:]) < 1.87 and \
		numpy.max(bestfit2[:]) < 1.925 and numpy.max(data2[:]+errors2[:]) < 1.925:
			ymax = 2.2
		elif numpy.max(bestfit1[:]) < 2.04 and numpy.max(data1[:]+errors1[:]) < 2.04 and \
		numpy.max(bestfit2[:]) < 2.1 and numpy.max(data2[:]+errors2[:]) < 2.1:
			ymax = 2.4
		elif numpy.max(bestfit1[:]) < 2.21 and numpy.max(data1[:]+errors1[:]) < 2.21 and \
		numpy.max(bestfit2[:]) < 2.275 and numpy.max(data2[:]+errors2[:]) < 2.275:
			ymax = 2.6
		elif numpy.max(bestfit1[:]) < 2.38 and numpy.max(data1[:]+errors1[:]) < 2.38 and \
		numpy.max(bestfit2[:]) < 2.45 and numpy.max(data2[:]+errors2[:]) < 2.45:
			ymax = 2.8
		elif numpy.max(bestfit1[:]) < 2.55 and numpy.max(data1[:]+errors1[:]) < 2.55 and \
		numpy.max(bestfit2[:]) < 2.625 and numpy.max(data2[:]+errors2[:]) < 2.625:
			ymax = 3.0
		elif numpy.max(bestfit1[:]) < 2.975 and numpy.max(data1[:]+errors1[:]) < 2.975 and \
		numpy.max(bestfit2[:]) < 3.0625 and numpy.max(data2[:]+errors2[:]) < 3.0625:
			ymax = 3.5
		elif numpy.max(bestfit1[:]) < 3.4 and numpy.max(data1[:]+errors1[:]) < 3.4 and \
		numpy.max(bestfit2[:]) < 3.5 and numpy.max(data2[:]+errors2[:]) < 3.5:
			ymax = 4.0
		elif numpy.max(bestfit1[:]) < 4.25 and numpy.max(data1[:]+errors1[:]) < 4.25 and \
		numpy.max(bestfit2[:]) < 4.375 and numpy.max(data2[:]+errors2[:]) < 4.375:
			ymax = 5.0
		elif numpy.max(bestfit1[:]) < 8.5 and numpy.max(data1[:]+errors1[:]) < 8.5 and \
		numpy.max(bestfit2[:]) < 8.75 and numpy.max(data2[:]+errors2[:]) < 8.75:
			ymax = 10.0
		elif numpy.max(bestfit1[:]) < 12.75 and numpy.max(data1[:]+errors1[:]) < 12.75 and \
		numpy.max(bestfit2[:]) < 13.125 and numpy.max(data2[:]+errors2[:]) < 13.125:
			ymax = 15.0
		else:
			ymax = 20.0
	if number_of_models==2:
		if numpy.max(bestfit1[:]) < 1.6 and numpy.max(data1[:]+errors1[:]) < 1.6 and \
		numpy.max(bestfit2[:]) < 1.75 and numpy.max(data2[:]+errors2[:]) < 1.75:
			ymax = 2.0
		elif numpy.max(bestfit1[:]) < 1.76 and numpy.max(data1[:]+errors1[:]) < 1.76 and \
		numpy.max(bestfit2[:]) < 1.925 and numpy.max(data2[:]+errors2[:]) < 1.925:
			ymax = 2.2
		elif numpy.max(bestfit1[:]) < 1.92 and numpy.max(data1[:]+errors1[:]) < 1.92 and \
		numpy.max(bestfit2[:]) < 2.1 and numpy.max(data2[:]+errors2[:]) < 2.1:
			ymax = 2.4
		elif numpy.max(bestfit1[:]) < 2.08 and numpy.max(data1[:]+errors1[:]) < 2.08 and \
		numpy.max(bestfit2[:]) < 2.275 and numpy.max(data2[:]+errors2[:]) < 2.275:
			ymax = 2.6
		elif numpy.max(bestfit1[:]) < 2.24 and numpy.max(data1[:]+errors1[:]) < 2.24 and \
		numpy.max(bestfit2[:]) < 2.45 and numpy.max(data2[:]+errors2[:]) < 2.45:
			ymax = 2.8
		elif numpy.max(bestfit1[:]) < 2.4 and numpy.max(data1[:]+errors1[:]) < 2.4 and \
		numpy.max(bestfit2[:]) < 2.625 and numpy.max(data2[:]+errors2[:]) < 2.625:
			ymax = 3.0
		elif numpy.max(bestfit1[:]) < 2.8 and numpy.max(data1[:]+errors1[:]) < 2.8 and \
		numpy.max(bestfit2[:]) < 3.0625 and numpy.max(data2[:]+errors2[:]) < 3.0625:
			ymax = 3.5
		elif numpy.max(bestfit1[:]) < 3.2 and numpy.max(data1[:]+errors1[:]) < 3.2 and \
		numpy.max(bestfit2[:]) < 3.5 and numpy.max(data2[:]+errors2[:]) < 3.5:
			ymax = 4.0
		elif numpy.max(bestfit1[:]) < 1.0 and numpy.max(data1[:]+errors1[:]) < 4.0 and \
		numpy.max(bestfit2[:]) < 4.375 and numpy.max(data2[:]+errors2[:]) < 4.375:
			ymax = 5.0
		elif numpy.max(bestfit1[:]) < 8.0 and numpy.max(data1[:]+errors1[:]) < 8.0 and \
		numpy.max(bestfit2[:]) < 8.75 and numpy.max(data2[:]+errors2[:]) < 8.75:
			ymax = 10.0
		elif numpy.max(bestfit1[:]) < 12.0 and numpy.max(data1[:]+errors1[:]) < 12.0 and \
		numpy.max(bestfit2[:]) < 13.125 and numpy.max(data2[:]+errors2[:]) < 13.125:
			ymax = 15.0
		else:
			ymax = 20.0
	if number_of_models==3:
		if numpy.max(bestfit1[:]) < 1.5 and numpy.max(data1[:]+errors1[:]) < 1.5 and \
		numpy.max(bestfit2[:]) < 1.75 and numpy.max(data2[:]+errors2[:]) < 1.75:
			ymax = 2.0
		elif numpy.max(bestfit1[:]) < 1.65 and numpy.max(data1[:]+errors1[:]) < 1.65 and \
		numpy.max(bestfit2[:]) < 1.925 and numpy.max(data2[:]+errors2[:]) < 1.925:
			ymax = 2.2
		elif numpy.max(bestfit1[:]) < 1.8 and numpy.max(data1[:]+errors1[:]) < 1.8 and \
		numpy.max(bestfit2[:]) < 2.1 and numpy.max(data2[:]+errors2[:]) < 2.1:
			ymax = 2.4
		elif numpy.max(bestfit1[:]) < 1.95 and numpy.max(data1[:]+errors1[:]) < 1.95 and \
		numpy.max(bestfit2[:]) < 2.275 and numpy.max(data2[:]+errors2[:]) < 2.275:
			ymax = 2.6
		elif numpy.max(bestfit1[:]) < 2.1 and numpy.max(data1[:]+errors1[:]) < 2.1 and \
		numpy.max(bestfit2[:]) < 2.45 and numpy.max(data2[:]+errors2[:]) < 2.45:
			ymax = 2.8
		elif numpy.max(bestfit1[:]) < 2.25 and numpy.max(data1[:]+errors1[:]) < 2.25 and \
		numpy.max(bestfit2[:]) < 2.625 and numpy.max(data2[:]+errors2[:]) < 2.625:
			ymax = 3.0
		elif numpy.max(bestfit1[:]) < 2.625 and numpy.max(data1[:]+errors1[:]) < 2.625 and \
		numpy.max(bestfit2[:]) < 3.0625 and numpy.max(data2[:]+errors2[:]) < 3.0625:
			ymax = 3.5
		elif numpy.max(bestfit1[:]) < 3.0 and numpy.max(data1[:]+errors1[:]) < 3.0 and \
		numpy.max(bestfit2[:]) < 3.5 and numpy.max(data2[:]+errors2[:]) < 3.5:
			ymax = 4.0
		elif numpy.max(bestfit1[:]) < 3.75 and numpy.max(data1[:]+errors1[:]) < 3.75 and \
		numpy.max(bestfit2[:]) < 4.375 and numpy.max(data2[:]+errors2[:]) < 4.375:
			ymax = 5.0
		elif numpy.max(bestfit1[:]) < 7.5 and numpy.max(data1[:]+errors1[:]) < 7.5 and \
		numpy.max(bestfit2[:]) < 8.75 and numpy.max(data2[:]+errors2[:]) < 8.75:
			ymax = 10.0
		elif numpy.max(bestfit1[:]) < 11.25 and numpy.max(data1[:]+errors1[:]) < 11.25 and \
		numpy.max(bestfit2[:]) < 13.125 and numpy.max(data2[:]+errors2[:]) < 13.125:
			ymax = 15.0
		else:
			ymax = 20.0
	if number_of_models==4:
		if numpy.max(bestfit1[:]) < 1.4 and numpy.max(data1[:]+errors1[:]) < 1.4 and \
		numpy.max(bestfit2[:]) < 1.75 and numpy.max(data2[:]+errors2[:]) < 1.75:
			ymax = 2.0
		elif numpy.max(bestfit1[:]) < 1.54 and numpy.max(data1[:]+errors1[:]) < 1.54 and \
		numpy.max(bestfit2[:]) < 1.925 and numpy.max(data2[:]+errors2[:]) < 1.925:
			ymax = 2.2
		elif numpy.max(bestfit1[:]) < 1.68 and numpy.max(data1[:]+errors1[:]) < 1.68 and \
		numpy.max(bestfit2[:]) < 2.1 and numpy.max(data2[:]+errors2[:]) < 2.1:
			ymax = 2.4
		elif numpy.max(bestfit1[:]) < 1.82 and numpy.max(data1[:]+errors1[:]) < 1.82 and \
		numpy.max(bestfit2[:]) < 2.275 and numpy.max(data2[:]+errors2[:]) < 2.275:
			ymax = 2.6
		elif numpy.max(bestfit1[:]) < 1.96 and numpy.max(data1[:]+errors1[:]) < 1.96 and \
		numpy.max(bestfit2[:]) < 2.45 and numpy.max(data2[:]+errors2[:]) < 2.45:
			ymax = 2.8
		elif numpy.max(bestfit1[:]) < 2.1 and numpy.max(data1[:]+errors1[:]) < 2.1 and \
		numpy.max(bestfit2[:]) < 2.625 and numpy.max(data2[:]+errors2[:]) < 2.625:
			ymax = 3.0
		elif numpy.max(bestfit1[:]) < 2.45 and numpy.max(data1[:]+errors1[:]) < 2.45 and \
		numpy.max(bestfit2[:]) < 3.0625 and numpy.max(data2[:]+errors2[:]) < 3.0625:
			ymax = 3.5
		elif numpy.max(bestfit1[:]) < 2.8 and numpy.max(data1[:]+errors1[:]) < 2.8 and \
		numpy.max(bestfit2[:]) < 3.5 and numpy.max(data2[:]+errors2[:]) < 3.5:
			ymax = 4.0
		elif numpy.max(bestfit1[:]) < 3.5 and numpy.max(data1[:]+errors1[:]) < 3.5 and \
		numpy.max(bestfit2[:]) < 4.375 and numpy.max(data2[:]+errors2[:]) < 4.375:
			ymax = 5.0
		elif numpy.max(bestfit1[:]) < 7.0 and numpy.max(data1[:]+errors1[:]) < 7.0 and \
		numpy.max(bestfit2[:]) < 8.75 and numpy.max(data2[:]+errors2[:]) < 8.75:
			ymax = 10.0
		elif numpy.max(bestfit1[:]) < 10.5 and numpy.max(data1[:]+errors1[:]) < 10.5 and \
		numpy.max(bestfit2[:]) < 13.125 and numpy.max(data2[:]+errors2[:]) < 13.125:
			ymax = 15.0
		else:
			ymax = 20.0
	return ymax


## Adds a "(SNIa)", "(SNcc)", or "(AGB)" after the name of each model (for the legend of the plot)
def rename_model(input_mod, list_models_SNIa, list_models_SNcc, list_models_AGB):
	if (input_mod in list_models_SNIa):
		input_mod = input_mod+" (SNIa)"
	elif (input_mod in list_models_SNcc):
		input_mod = input_mod+" (SNcc)"
	else:
		input_mod = input_mod+" (AGB)"
	return input_mod
	
	
## Labels the SN fraction accordingly (for the legend of the plot)
def label_SNfrac(input_mod1, input_mod2, list_models_SNIa, list_models_SNcc, list_models_AGB):
	if (input_mod1 in list_models_SNIa) and (input_mod2 in list_models_SNcc):
		fraction_label = "$\\frac{\mathrm{SNIa}}{\mathrm{SNIa+SNcc}} = $"
	elif (input_mod1 in list_models_SNcc) and (input_mod2 in list_models_SNIa):
		fraction_label = "$\\frac{\mathrm{SNcc}}{\mathrm{SNIa+SNcc}} = $"
	else:
		fraction_label = "$\\frac{\mathrm{model1}}{\mathrm{all\_models}} = $"
	return fraction_label
	


## Selects a color for models and avoid conflicts in case of several models from the same type  (for the legend of the plot)
def set_colors_1(input_mod_1, default_color_SNIa, default_color_SNcc, default_color_AGB):
	list_models_SNIa = numpy.loadtxt("list_models_SNIa.txt", dtype="str").tolist()
	list_models_SNcc = numpy.loadtxt("list_models_SNcc.txt", dtype="str").tolist()
	list_models_AGB = numpy.loadtxt("list_models_AGB.txt", dtype="str").tolist()

	if (input_mod_1 in list_models_SNIa):
		color1 = default_color_SNIa
	elif (input_mod_1 in list_models_SNcc):
		color1 = default_color_SNcc
	else:
		color1 = default_color_AGB
	return color1


## Selects a color for models and avoid conflicts in case of several models from the same type  (for the legend of the plot)
def set_colors_2(input_mod_1, input_mod_2, default_color_SNIa, default_color_SNcc, default_color_AGB):
	list_models_SNIa = numpy.loadtxt("list_models_SNIa.txt", dtype="str").tolist()
	list_models_SNcc = numpy.loadtxt("list_models_SNcc.txt", dtype="str").tolist()
	list_models_AGB = numpy.loadtxt("list_models_AGB.txt", dtype="str").tolist()

	color2 = ["aaa", "bbb"]
	color2[0] = set_colors_1(input_mod_1, default_color_SNIa, default_color_SNcc, default_color_AGB)

	if (input_mod_2 in list_models_SNIa):
		if (input_mod_1 in list_models_SNIa):
			color2[1] = "darkblue"
		else:
			color2[1] = default_color_SNIa
	elif (input_mod_2 in list_models_SNcc):
		if (input_mod_1 in list_models_SNcc):
			color2[1] = "darkviolet"
		else:
			color2[1] = default_color_SNcc
	elif (input_mod_2 in list_models_AGB):
		if (input_mod_1 in list_models_AGB):
			color2[1] = "orange"
		else:
			color2[1] = default_color_AGB
	return color2


## Selects a color for models and avoid conflicts in case of several models from the same type  (for the legend of the plot)
def set_colors_3(input_mod_1, input_mod_2, input_mod_3, default_color_SNIa, default_color_SNcc, default_color_AGB):
	list_models_SNIa = numpy.loadtxt("list_models_SNIa.txt", dtype="str").tolist()
	list_models_SNcc = numpy.loadtxt("list_models_SNcc.txt", dtype="str").tolist()
	list_models_AGB = numpy.loadtxt("list_models_AGB.txt", dtype="str").tolist()

	color3 = ["aaa", "bbb", "ccc"]
	color3[0:2] = set_colors_2(input_mod_1, input_mod_2, default_color_SNIa, default_color_SNcc, default_color_AGB)

	if (input_mod_3 in list_models_SNIa):
		if (input_mod_1 in list_models_SNIa) or (input_mod_2 in list_models_SNIa):
			color3[2] = "green"
		else:
			color3[2] = default_color_SNIa
	elif (input_mod_3 in list_models_SNcc):

		if (input_mod_1 in list_models_SNcc) or (input_mod_2 in list_models_SNcc):
			color3[2] = "mediumpurple"
		else:
			color3[2] = default_color_SNcc
	else:
		if (input_mod_1 in list_models_AGB) or (input_mod_2 in list_models_AGB):
			color3[2] = "red"
		else:
			color3[2] = default_color_AGB
	return color3


## Selects a color for models and avoid conflicts in case of several models from the same type  (for the legend of the plot)
def set_colors_4(input_mod_1, input_mod_2, input_mod_3, input_mod_4):
	list_models_SNIa = numpy.loadtxt("list_models_SNIa.txt", dtype="str").tolist()
	list_models_SNcc = numpy.loadtxt("list_models_SNcc.txt", dtype="str").tolist()
	list_models_AGB = numpy.loadtxt("list_models_AGB.txt", dtype="str").tolist()

	color4 = ["aaa", "bbb", "ccc", "ddd"]
	color4[0:3] = set_colors_3(input_mod_1, input_mod_2, input_mod_3, default_color_SNIa, default_color_SNcc, default_color_AGB)

	if (input_mod_4 in list_models_SNIa):
		if (input_mod_1 in list_models_SNIa) or (input_mod_2 in list_models_SNIa) \
		or (input_mod_3 in list_models_SNIa):
			color4[3] = "gray"
		else:
			color4[3] = default_color_SNIa
	elif (input_mod_3 in list_models_SNcc):
		if (input_mod_1 in list_models_SNcc) or (input_mod_2 in list_models_SNcc) \
		or (input_mod_3 in list_models_SNcc):
			color4[3] = "gray"
		else:
			color4[3] = default_color_SNcc
	else:
		if (input_mod_1 in list_models_AGB) or (input_mod_2 in list_models_AGB) \
		or (input_mod_3 in list_models_AGB):
			color4[3] = "gray"
		else:
			color4[3] = default_color_AGB
	return color4


   
main()
