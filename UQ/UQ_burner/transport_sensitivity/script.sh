#for file in stagnation_flame_uiuc20sp_m25.00_phi1.00_S????_coeff0.50.csv
#do
#	ireact=$(head -n 1 $file)
#	react=$(head -n 4 $file | tail -n 1)
#	flux=$(head -n 3 $file | tail -n 1)
#	temp=$(head -n 2 $file | tail -n 1)
#	echo $flux $temp $ireact $react
#done

rm -rf output
mkdir -p ./output
mech="wang99_51sp"
phi_array=("0.55" "0.70" "0.85" "1.00" "1.15" "1.30") 
param_array=("diameter" "well-depth" "dipole" "polarizability" "rotational-relaxation")
species_array=("C2H4" "H2" "H" "O2" "O" "H2O" "CO" "CO2" "OH" "HCO" "HO2" "H2O2" "C2H3" "C2H2" "CH4" "CH3" "CH2" "CH2O" "CH2CHO" "N2")
coeff_array=("0.90" "1.10")
for coeff in "${coeff_array[@]}"
do
    for parameter in "${param_array[@]}"
    do
        for species in "${species_array[@]}"
        do
            for phi in "${phi_array[@]}"
            do
                file=./flux/stagnation_flame_${mech}_m25.00_phi${phi}_${parameter}${species}_coeff${coeff}.dat 
                if [ -f "$file" ] 
        	    then
    	            spc=$(head -n 1 $file)
    	            flux=$(head -n 3 $file | tail -n 1)
    	            temp=$(head -n 2 $file | tail -n 1)
                    old_param=$(head -n 4 $file | tail -n 1)
                    new_param=$(head -n 5 $file | tail -n 1)
                    parameter=$(head -n 6 $file | tail -n 1)
    	            echo $phi $flux $temp $old_param $new_param \"$parameter\" >> ./output/${parameter}${species}_${coeff}.dat
                fi 
            done
        done
    done
done
