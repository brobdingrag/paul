
# Install ped-sim
git clone https://github.com/williamslab/ped-sim.git

cd ped-sim

# Make
cp Makefile-gsl Makefile
make

# Download European genetic map
wget https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/raw/master/Refined_EUR_genetic_map_b37.tar.gz
tar xvzf Refined_EUR_genetic_map_b37.tar.gz
printf "#chr\tpos\tmale_cM\tfemale_cM\n" > refined_mf.simmap
for chr in {1..22}; do
  paste Refined_EUR_genetic_map_b37/male_chr$chr.txt Refined_EUR_genetic_map_b37/female_chr$chr.txt \
    | awk -v OFS="\t" 'NR > 1 && $2 == $6 {print $1,$2,$4,$8}' \
    | sed 's/^chr//' >> refined_mf.simmap;
done

mv refined_mf.simmap refined_mf_X.simmap
awk 'NR > 1 { print $1,$2,"0.0",$4 }' Refined_EUR_genetic_map_b37/female_chrX.txt \
    | sed 's/^chr//' >> refined_mf_X.simmap


# Define the filename
DEF_FILENAME="father-child.def"

# Create the file with the specified contents using a here-document
cat << EOF > "$DEF_FILENAME"
# Define pedigree name, number of copies, and number of generations
def father-child 1 2 M

# Generation 1: Print the father
1 1

# Generation 2: Print children in 1 branch
2 10000 1
EOF

echo "Running ped-sim..."

# Run the simulation with interference and showing breakpoints
./ped-sim -d father-child.def -m refined_mf_X.simmap -o father_child_bp --intf interfere/nu_p_campbell_X.tsv --bp --seed 42

# Get the crossover positions
python3 save_meioses.py

