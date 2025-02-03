rm -rf output
mkdir -p ./output
phi_array=("0.55" "0.70" "0.85" "1.00" "1.15" "1.30")
coeff_array=("0.50" "2.00")
for coeff in "${coeff_array[@]}"
do
    for rxn in {0000..0077}
    do
        for phi in "${phi_array[@]}"
        do
            file=./flux/stagnation_flame_uiuc20sp_m25.00_phi${phi}_S${rxn}_coeff${coeff}.csv
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
