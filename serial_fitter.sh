#! /bin/bash

############ SERIAL FITTER ############
#######################################

# Important: Please set "plotting_mode = False" in abunfit.py before running this script...

#(chi^2) given here is actually (chi^2)*100000 !!


echo "10000000000000" > old_chi.txt
rm summary_models.txt
rm summary_models_sorted.txt




for SNcc in No13_SNcc_0 No13_SNcc_0.001 No13_SNcc_0.004 No13_SNcc_0.008 No13_SNcc_0.02 Su16_N20 Su16_W18
do



for SNIa in Se13_N1 Se13_N3 Se13_N5 Se13_N10 Se13_N20 Se13_N40 Se13_N100H Se13_N100 Se13_N100L Se13_N150 Se13_N200 Se13_N300C Se13_N1600 Se13_N1600C Se13_N100_Z0.5 Se13_N100_Z0.1 Se13_N100_Z0.01 Fi14_N1def Fi14_N3def Fi14_N5def Fi14_N10def Fi14_N20def Fi14_N40def Fi14_N100Hdef Fi14_N100def Fi14_N100Ldef Fi14_N150def Fi14_N200def Fi14_N300Cdef Fi14_N1600def Fi14_N1600Cdef Le18_050-1-c3-1P Le18_100-1-c3-1P Le18_100-0-c3 Le18_100-0.1-c3 Le18_100-0.5-c3 Le18_100-1-c3 Le18_100-2-c3 Le18_100-3-c3 Le18_100-5-c3 Le18_300-1-c3-1P Le18_300-0-c3 Le18_300-0.1-c3 Le18_300-0.5-c3 Le18_300-1-c3 Le18_300-2-c3 Le18_300-3-c3 Le18_300-5-c3 Le18_500-1-c3-1P Le18_500-0-c3 Le18_500-0.1-c3 Le18_500-0.5-c3 Le18_500-1-c3 Le18_500-2-c3 Le18_500-3-c3 Le18_500-5-c3
do



for alpha in -2.35
#for alpha in 0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2 -1.4 -1.6 -1.8 -2.0 -2.1 -2.20 -2.25 -2.30 -2.35 -2.40 -2.45 -2.50 -2.6 -2.7 -2.8 -2.9 -3.0
do



echo "===================="
echo ${SNcc} ${SNIa}
./abunfit.py data/CHEERS_Mernier18b.out 2 ${alpha} ${SNcc} ${SNIa}


old_chi=$( sed -n '1p' old_chi.txt )
chi=$( sed -n '1p' currentchi.txt )

echo ${chi} ${alpha} ${SNcc} ${SNIa} >> summary_models.txt

if [ ${chi} -lt ${old_chi} ]; then
echo ${chi} > old_chi.txt
echo ${SNIa} > best_model.txt
echo ${SNcc} >> best_model.txt
echo ${alpha} >> best_model.txt
echo ${chi} >> best_model.txt
fi



done
done
done




SNIa=$( sed -n '1p' best_model.txt )
SNcc=$( sed -n '2p' best_model.txt )
alpha=$( sed -n '3p' best_model.txt )

./abunfit.py data/CHEERS_Mernier18b.out 2 ${alpha} ${SNcc} ${SNIa}


./sort_summary_models.py summary_models.txt
mv summary_models_sorted.txt summary_models.txt

