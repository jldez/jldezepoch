#!/bin/bash

# shoot laser on stuff
material='silica'
lambda=800 	#nm
fluence=2.0 	#J/cm2
tau=10 		#fs
N=1

# Simulation parameters
depth=120			#nm
width_para=80		#nm
width_perp=80		#nm
res_depth=4 		#nm
res_width=20 		#nm
nb_snapshots=100

# Create directories and files from the generic
dir=$material"_"$lambda"_"$fluence"_"$tau
mkdir $dir
cp -r generic/* $dir
cp materialParams/$material.sh $dir

# Import material parameters
source materialParams/$material.sh

# Define parameters for ionization rate and calculate
lambda_line_keldysh="lambd = "$lambda".0"
bandGap_line_keldysh="bandGap_eV = "$bandGap
mReduced=$(echo "scale=6 ; 1.0 / $effective_mass_hole + 1.0 / $effective_mass_hole" | bc)
mReduced=$(echo "scale=6 ; 1.0 / $mReduced" | bc)
mReduced_line_keldysh="mReduced = "$mReduced""
densiteAtome_line_keldysh="densiteAtome = "$atom_density
perl -i -pe "s/.*/$lambda_line_keldysh/ if $. == 36" $dir"/keldyshRate.py"
perl -i -pe "s/.*/$bandGap_line_keldysh/ if $. == 45" $dir"/keldyshRate.py"
perl -i -pe "s/.*/$mReduced_line_keldysh/ if $. == 54" $dir"/keldyshRate.py"
perl -i -pe "s/.*/$densiteAtome_line_keldysh/ if $. == 57" $dir"/keldyshRate.py"
# python $dir"/keldyshRate.py" #uncomment to generate ionization table

# Produce the input deck for N=1
nx_line_deck="	nx = $(($(($depth+50))/$res_depth))"
ny_line_deck="	ny = $(($width_perp/$res_width))"
nz_line_deck="	nz = $(($width_para/$res_width))"
t_end_line_deck="	t_end = "$tau" * 3.001 * femto"
x_max_line_deck="	x_max = $(($depth+50))e-9"
y_max_line_deck="	y_max = "$width_perp"e-9"
z_max_line_deck="	z_max = "$width_para"e-9"
path_line_deck="	path_ionRate_table = "$dir
intensity_line_deck="	intensity_w_cm2 = "$fluence" \/"$tau" * 1.595769e15"
lambda_line_deck="	lambda = "$lambda"e-9"
t_profile_line_deck="	t_profile = gauss(time,"$tau"*1.5*femto,"$tau"*femto)"
path_mm_line_deck="	path_medium_mask = "$dir"\/pulse1"
path_ri_line_deck="	path_rho_incubation = "$dir"\/pulse1"
resonance_line_deck="	resonance = "$resonance
gamma_line_deck="	gamma = "$gamma
chi1_line_deck="	chi1 = "$chi1
chi2_line_deck="	chi2 = "$chi2
chi3_line_deck="	chi3 = "$chi3
gamma_D_line_deck="	gamma_D = "$gamma_D
bandGap_line_deck="	bandGap_drude = "$bandGap"*ev"
avalanche_line_deck="	avalanche_factor = "$avalanche_factor
atom_density_line_deck="	atom_density = "$atom_density
atom_cross_section_line_deck="	atom_cross_section = "$atom_cross_section
effective_mass_electron_line_deck="	effective_mass_electron = "$effective_mass_electron
effective_mass_hole_line_deck="	effective_mass_hole = "$effective_mass_hole
lambda_to_omega=1.88365e18
omega_laser_line_deck="	omega_laser = "$lambda_to_omega" \/ "$lambda
amp_laser_line_deck="	amp_laser = sqrt("$fluence"\/"$tau") * 6.5199274e10"
recombination_rate_line_deck="	recombination_rate = "$recombination_rate
snap_time_line_deck="	dt_snapshot = "$tau" * 3 \/ "$nb_snapshots" * femto"

perl -i -pe "s/.*/$nx_line_deck/ if $. == 2" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$ny_line_deck/ if $. == 3" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$nz_line_deck/ if $. == 4" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$t_end_line_deck/ if $. == 7" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$x_max_line_deck/ if $. == 13" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$y_max_line_deck/ if $. == 15" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$z_max_line_deck/ if $. == 17" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$path_line_deck/ if $. == 23" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$intensity_line_deck/ if $. == 41" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$lambda_line_deck/ if $. == 42" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$t_profile_line_deck/ if $. == 43" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$path_mm_line_deck/ if $. == 47" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$path_ri_line_deck/ if $. == 48" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$resonance_line_deck/ if $. == 58" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$gamma_line_deck/ if $. == 59" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$chi1_line_deck/ if $. == 60" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$chi2_line_deck/ if $. == 61" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$chi3_line_deck/ if $. == 62" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$gamma_D_line_deck/ if $. == 63" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$bandGap_line_deck/ if $. == 64" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$avalanche_line_deck/ if $. == 65" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$atom_density_line_deck/ if $. == 66" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$atom_cross_section_line_deck/ if $. == 67" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$effective_mass_electron_line_deck/ if $. == 68" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$effective_mass_hole_line_deck/ if $. == 69" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$recombination_rate_line_deck/ if $. == 70" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$omega_laser_line_deck/ if $. == 71" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$amp_laser_line_deck/ if $. == 72" $dir"/pulse1/input.deck"
perl -i -pe "s/.*/$snap_time_line_deck/ if $. == 76" $dir"/pulse1/input.deck"

# Input decks for N>1
for (( n=2; n<=N+1; n++ ))
do
	mkdir $dir"/pulse"$n 
	cp $dir"/pulse$(($n-1))/input.deck" $dir"/pulse"$n"/input.deck"
	path_mm_line_deck="	path_medium_mask = "$dir"\/pulse"$n
	path_ri_line_deck="	path_rho_incubation = "$dir"\/pulse"$n
	perl -i -pe "s/.*/$path_mm_line_deck/ if $. == 47" $dir"/pulse"$n"/input.deck"
	perl -i -pe "s/.*/$path_ri_line_deck/ if $. == 48" $dir"/pulse"$n"/input.deck"
done

# Produce interpulse.py
trap_density_line_interpulse="trap_density = "$trap_density
perl -i -pe "s/.*/$trap_density_line_interpulse/ if $. == 16" $dir"/interpulse.py"
xi_line_interpulse="xi = "$xi
perl -i -pe "s/.*/$xi_line_interpulse/ if $. == 17" $dir"/interpulse.py"
rho_sat_line_interpulse="rho_sat = "$atom_density
perl -i -pe "s/.*/$rho_sat_line_interpulse/ if $. == 18" $dir"/interpulse.py"
Eb_line_interpulse="Eb = "$Eb
perl -i -pe "s/.*/$Eb_line_interpulse/ if $. == 19" $dir"/interpulse.py"
bandGap_line_interpulse="bandGap = "$bandGap
perl -i -pe "s/.*/$bandGap_line_interpulse/ if $. == 20" $dir"/interpulse.py"

# Run
for (( n=1; n<=N; n++ ))
do
	mpirun -np 1 ./bin/epoch3d <<< $dir"/pulse"$n
	# python $dir/interpulse.py $n $nb_snapshots $dir
done