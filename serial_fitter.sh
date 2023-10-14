#! /bin/bash

############ SERIAL FITTER ############
#######################################


#(chi^2) given here is actually (chi^2)*100000 !!


echo "10000000000000" > old_chi.txt
rm fitting_results_sorted.txt



for alpha in -2.35
do
 
	for SNcc in No13_SNcc_0 No13_SNcc_0.001 No13_SNcc_0.004 No13_SNcc_0.008 No13_SNcc_0.02 Su16_N20 Su16_W18
	do

		for SNIa in Se13_N1 Se13_N3 Se13_N5 Se13_N10 Se13_N20 Se13_N40 Se13_N100H Se13_N100 Se13_N100L Se13_N150 Se13_N200 Se13_N300C Se13_N1600 Se13_N1600C Se13_N100_Z0.5 Se13_N100_Z0.1 Se13_N100_Z0.01
		do


			echo "===================="
			echo ${SNcc} ${SNIa}
			python abunfit.py data/CHEERS_Mernier18b.dat 2 ${alpha} ${SNcc} ${SNIa} --disable-plot


			old_chi=$( sed -n '1p' old_chi.txt )
			chi=$( sed -n '1p' currentchi.txt )

			echo ${chi} ${alpha} ${SNIa} ${SNcc} >> fitting_results.txt


		done
	done
done



python sort_fitting_results.py fitting_results.txt fitting_results_sorted.txt
rm fitting_results.txt
