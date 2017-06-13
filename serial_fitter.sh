#! /bin/bash
#chi given here is actual chi^2*100000 !!


echo "10000000000000" > old_chi.txt
rm summary_models.txt
rm summary_models_sorted.txt




for SNcc in N13_SNcc_0 N13_SNcc_0.001 N13_SNcc_0.004 N13_SNcc_0.008 N13_SNcc_0.02

do




for SNIa in Seitenzahl_N1 Seitenzahl_N3 Seitenzahl_N5 Seitenzahl_N10 Seitenzahl_N20 Seitenzahl_N40 Seitenzahl_N100H Seitenzahl_N100 Seitenzahl_N100L Seitenzahl_N150 Seitenzahl_N200 Seitenzahl_N300C Seitenzahl_N1600 Seitenzahl_N1600C Seitenzahl_N100_Z0.5 Seitenzahl_N100_Z0.1 Seitenzahl_N100_Z0.01 Fink_N1def Fink_N3def Fink_N5def Fink_N10def Fink_N20def Fink_N40def Fink_N100Hdef Fink_N100def Fink_N100Ldef Fink_N150def Fink_N200def Fink_N300Cdef Fink_N1600def Fink_N1600Cdef

do



for alpha in -2.35
#for alpha in 0.0 -0.2 -0.4 -0.6 -0.8 -1.0 -1.2 -1.4 -1.6 -1.8 -2.0 -2.1 -2.20 -2.25 -2.30 -2.35 -2.40 -2.45 -2.50 -2.6 -2.7 -2.8 -2.9 -3.0
do



echo "===================="
echo ${SNcc} ${SNIa}
./abunfit.py cluster/average_CHEERS.out 2 ${alpha} ${SNcc} ${SNIa}


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

./abunfit.py cluster/average_CHEERS.out 2 ${alpha} ${SNcc} ${SNIa}


./sort_summary_models.py summary_models.txt
mv summary_models_sorted.txt summary_models.txt

