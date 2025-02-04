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
phi_array=("0.55" "0.70" "0.85" "1.00" "1.15" "1.30") 
coeff_array=("0.90" "1.10")
mech="wang99_51sp"
for coeff in "${coeff_array[@]}"
do
    for rxn in {0000..0026}
    do
        for phi in "${phi_array[@]}"
        do
            file=./flux/stagnation_flame_${mech}_m25.00_phi${phi}_S${rxn}_coeff${coeff}.csv 
            if [ -f "$file" ] 
    	    then
	        ireact=$(head -n 1 $file)
	        react=$(head -n 4 $file | tail -n 1)
	        flux=$(head -n 3 $file | tail -n 1)
	        temp=$(head -n 2 $file | tail -n 1)
	        echo $phi $flux $temp $ireact \"$react\" >> ./output/reaction_${rxn}_${coeff}.dat
            fi 
        done
    done
done
