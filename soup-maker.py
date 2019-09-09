#!/usr/bin/env python3

import argparse, os, glob, math
import subprocess
import numpy as np

export GMX_MAXBACKUP=-1

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmx", "--gromacs_command",  type=str,  default="gmx")
    # parser.add_argument("-cr_n", "--n", help="crowder index",  type=int,  default=0)
    
    args = parser.parse_args()
    return args

box_size_initial=35		# in nm
box_size_step=2	# nm

droplet_thickness_initial=0.6

droplet_shell_precision=2	# percentage as integer

crowder_total=13
fraction=0.30
ionic_strength=0.150	# M

make_droplets=0
make_box=1
make_top=1
make_neutral=1
make_is=1
check_fraction=1

path=os.getcwd()
path_crowders=path
path_mdp='mdp'

mdp_min='min.mdp'

path_output='soup_4'

# number of gro files for each crowder
num_crowder[1]=6
num_crowder[2]=7
num_crowder[3]=2
num_crowder[4]=1
num_crowder[5]=3
num_crowder[6]=1
num_crowder[7]=1
num_crowder[8]=1
num_crowder[9]=5
num_crowder[10]=1436
num_crowder[11]=144
num_crowder[12]=225
num_crowder[13]=255
crowder_type_arr = []
crowder_mass = []
# number of gro files for each crowder


print('Preliminary steps\n')
print('\tGathering data about each crowder\n')

crowder=1
for crowder in range(1, crowder_total+1):
	# Creates a temporary index file to check what type of molecule is the crowder
	print( 'q\n | gmx make_ndx -f crowder_n%s/npt.1.gro -o temporary.ndx >/dev/null 2>/dev/null ' %(crowder))
	with open('q.txt', 'w') as outfile:
            outfile.write("q\n")
    os.system('%s make_ndx -f crowder_n%s/npt.1.gro -o temporary.ndx < q.txt  >/dev/null 2>/dev/null ' % (gromacs_command, str(n), str(n) ))
    os.remove('q.txt')
    with open('temporary.ndx', 'r') as infile:
        l=[]
	# Is this crowder a protein, RNA or metabolite?
        for line in infile:
            if line.find('[') > -1:
                l.append(line.strip())
        if '[ Protein ]' in l:
            crowder_type = "Protein"
        elif '[ RNA ]' in l:
            crowder_type = "RNA"
        elif '[ Other ]' in l:
            crowder_type = "Other"            
        else:
            print("Unknown type for crowder_n%s" % str(n))

	crowder_type_arr.append(crowder_type)
	
	# calculate the mass of each crowder (considering all chains) from the merged itp files
    with open('crowder_n%s/crowder_n%s.itp' % (str(n), str(n)), 'r') as infile:
        mass = []
        for line in infile:
            if line.find('qtot'):
                mass.append(float(line[7]))
    crowder_mass.append(sum(mass))

print('\tSorting the crowders by mass')
indices_sorted = np.flip(np.argsort(crowder_mass))+1

print('\tCalculating the biomolecular mass of the soup')
total_biomol_mass=0
for crowder in range(1, crowder_total+1):
    total_biomol_mass = total_biomol_mass + num_crowder[crowder]*crowder_mass[crowder-1]


print('\t\tTotal biomolecular mass =%0.3f' % (total_nonwater_mass))

# echo $fraction | awk '{printf "\tCalculating the number of waters necessary to achieve the biomolecular fraction of %.0f %%\n", $1*100}'

# water_necessary=$(echo $total_biomol_mass $fraction | awk '{printf "%.0f\n", (($1/$2)-$1)/18 }')

# total_water_mass=$(echo $water_necessary | awk '{printf "%.2f\n",$1*18}')

# echo -e "\t\t"'It is necessary to add '"$water_necessary"' water molecules to the final box'
# echo -e "\t\t"'That corresponds to '"$total_water_mass"' Da'
# echo
# echo
# echo



if make_droplets == 1:
    print('Making droplets for each crowder')
    print('\tFinding the optimum droplet shell thickness which comprises a total of %s water molecules for all crowders\n' % str(water_necessary))
    print('\tThe maximum difference that is allowed between the number of water molecules inside the droplets and the number of molecules we need is %s %%' % str(droplet_shell_precision))
# if [ $make_droplets -eq 1 ]; then
	
    ####Zaza: convert the line below
# 	echo 'thickness[nm] biomol-fract[%] added-perc[%] error-perc[%]' | awk '{printf "\t\t%16s\t%16s\t%16s\t%16s\n",$1,$2,$3,$4}'
	
	droplet_thickness_step=0.1	# nm
	droplet_thickness=droplet_thickness_initial		# nm
	droplet_thickness_prev=droplet_thickness
	
	ok=0
    while ok  == 0 :
        water_total_added = 0
        for crowder in range(1, crowder_total+1):
            os.chdir('crowder_n%s' % (str(crowder)))
			# Generate a box with PBC corrections
            os.system('echo System | %s trjconv -f npt.1.gro -s npt.1.tpr -pbc atom -ur compact -o eq_pbc.gro >/dev/null 2>/dev/null' % (args.gromacs_command))
          # Create the index_eq.ndx file containing Crowder + Waters+Ions within a shell around the crowder
            os.system('%s select -f eq_pbc.gro -on index_droplet.ndx -select "(group %s ) or (same residue as (group "SOL" and within %s of group %s) ) or ( group "Ion" and within %s of group %s )" -s npt.1.tpr >/dev/null 2>dev/null ' % (args.gromacs_command, crowder_type[crowder], str(droplet_thickness), crowder_type[crowder], str(droplet_thickness), crowder_type[crowder]))
          # Create the gro file with the current droplet          
            os.system('%s trjconv -f eq_pbc.gro -s npt.1.tpr -o temporary.gro -n index_droplet.ndx >/dev/null 2>dev/null')
            with open('droplet.gro', 'w') as outfile:
                outfile.write('droplet, crowder_n%s' % str(crowder))

            with open('temporary.gro', 'r') as infile:
                remove_k=0
                remove_cl=0
                for i, line in enumerate(infile):
                    if i == 1:
                        atoms = int(line[0])
                    if line.find('K      K'):
                        remove_k = remove_k+1
                    if line.find('CL      CL'):
                        remove_cl = remove_cl+1
            with open('droplet.gro', 'a') as outfile:
                outfile.write('%s' % str(atoms-remove_cl-remove_k))

			with open('droplet.gro', 'a') as outfile:
                with open('temporary.gro', 'r') as infile:
                    for i, line in enumerate(infile):
                        if i > 1:
                            if line.find('K        K') < 0 or line.find('CL        CL') < 0:
                                outfile.write(line)
			
            os.remove('index_droplet.ndx')
			# Number of waters in the droplet
            water_added_by_one_droplet=0
            with open('droplet.gro', 'r') as infile:
                for line in infile:
                    if line.find(" OW"):
                        water_added_by_one_droplet=water_added_by_one_droplet+1
			
			# Total waters to be added by this crowder
			water_added_by_this_crowder=water_added_by_one_droplet*num_crowder[crowder]
			
			# Update the current number of waters in the soup
			water_total_added=water_total_added+water_added_by_this_crowder
			os.chdir('..')

        water_added_perc = (float(water_total_added)/float(water_necessary)) * 100
        error_perc = water_added_perc - 100
		error_perc_abs = -1*error_perc

		# Fraction 
        biomol_fraction = float(total_biomol_mass) / (float(total_biomol_mass) + float(water_total_added)*18.)
		
		# Output
        print("\t\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n" % (float(droplet_thickness), float(biomol_fraction)*100., float(water_added_perc), float(error_perc)))		
		
        if float(error_perc_abs) <= float(droplet_shell_precision) :
            ok=1
        else:
          # If we have missing waters, increase the shell size and try again
            if float(water_added_perc < 100) :
                droplet_thickness_prev = float(droplet_thickness)
            else:
                droplet_thickness=float(droplet_thickness_prev)
                droplet_thickness_step=float(droplet_thickness_prev)/10.
            droplet_thickness=float(droplet_thickness) + float(droplet_thickness_step)


    print('\tDroplet shell thickness converged to %s nm.' % (str(droplet_thickness) ))
    print('\t %d water molecules will be added to the box. This is %0.2f \% of the precise number required, which is within the requested error margin of %0.2f \%' % (water_total_added, water_added_perc, droplet_shell_precision) )
    print('\n')     
    print('\n')     
    print('\n')     

    os.chdir(path)


if make_box == 1:
    print('\tFinding the smallest cubic box in which the soup fits')
    # Sorts the list of crowders according to their number of atoms from biggest to smallest
    # This is done because its easier to fit all of them in a box if we start with the biggest
    os.mkdir(path_output)
    # escape condition
    ok=0
    box_size = float(box_size_initial)
    while True:
        # we start with an empty cube with sides $box_size nm
        with open('%s/box.gro' % (path_output), 'w') as outfile:
            outfile.write('CYTOPLASM MODEL\n')
            outfile.write('0\n')
            outfile.write('%s\t%s\t%s' % (str(box_size), str(box_size), str(box_size)))
        print('\t\tSide Lenght = %s nm' % (str(box_size)))

        # for each crowder ...
        crowder_type_count = 1
        for crowder in indices_sorted:
            print('\t\t\tTrying to add %s copies of crowder number %s ... ' % (num_crowder[crowder], str(crowder)))
            os.system('%s insert-molecules -f %s/box.gro -ci crowder_n%s/droplet.gro -o %s/box.gro -nmol %s -try 10 >out 2>err' % (path_output, str(crowder), path_output, str(num_crowder[crowder])))
            with open('err', 'r') as infile:
                for line in infile:
                    if line.find('Added') > -1:
                        number_added = int(line.strip()[1])
            os.remove("err")
            os.remove("out")
            print('added %s ...' % (str(number_added)))
            # If we added as many as we wanted ...
            if number_added == num_crowder[crowder]:
                print('ok')
                if crowder_type_count == crowder_total:
                    # Stops if this is the last crowder
                    ok=1
                    break
                else:
                    crowder_type_count = crowder_type_count +1
            else:
                # If we couln't add as many as we wanted...
                box_size = box_size + box_size_step
                print('didn't fit)
                # Get out of the loop and start again
                break
        #If there are no more crowders and everything was ok, stop
        if ok == 1:
            break

    print('\tWriting the .gro file for the soup ... ')
    # Splits the current box by its constituents 
    K_temp, CL_temp, MG_temp, SOL_temp = [],[],[],[]
    k_counter = 0
    cl_counter = 0
    mg_counter = 0
    sol_counter = 0
    with open('%s/box.gro' % (path_output), 'r') as infile:
        with open('biomol_temp', 'w') as outfile:
            for line in infile:
                if line.find("K        K") > -1 :
                    K_temp.append(line.strip())
                    k_counter = k_counter + 1
                if line.find("CL        CL") > -1 :
                    CL_temp.append(line.strip())
                    cl_counter = cl_counter + 1
                if line.find("MG        MG") > -1 :
                    MG_temp .append(line.strip())
                    mg_counter = mg_counter + 1
                if line.find("SOL ") > -1 :
                    SOL_temp.append(line.strip())
                    sol_counter = sol_counter + 1
                else:
    # 	cat "$path_output"/box.gro | grep -v "K        K" | grep -v "CL      CL" | grep -v "MG      MG" | grep -v "SOL " | tail -n+3 | head -n-1 > biomol_temp
                    outfile.write(line)
    
    with open('k_temp', 'w') as outfile:
        for line in K_temp:
            outfile.write(line)
    with open('cl_temp', 'w') as outfile:
        for line in CL_temp:
            outfile.write(line)
    with open('mg_temp', 'w') as outfile:
        for line in MG_temp:
            outfile.write(line)
    with open('sol_temp', 'w') as outfile:
        for line in SOL_temp:
            outfile.write(line)
    
	
# 	number_total=$(cat biomol_temp SOL_temp K_temp CL_temp MG_temp | wc -l | awk '{print $1}')

	# Writes a new gro file in the correct order
    with open('%s/box_ordered.gro' % (path_output), 'r') as infile:
        infile.write("SOUP\n")
        infile.write("%s" % str(number_total))
        with open('biomol_temp', 'r') as outfile:
            for line in outfile:
                infile.write(line)
        with open('k_temp', 'r') as outfile:
            for line in outfile:
                infile.write(line)
        with open('cl_temp', 'r') as outfile:
            for line in outfile:
                infile.write(line)
        with open('mg_temp', 'r') as outfile:
            for line in outfile:
                infile.write(line)
        with open('sol_temp', 'r') as outfile:
            for line in outfile:
                infile.write(line)
        with open('%s/box.gro' % (path_output), 'r') as outfile:
            for line in outfile:
                True
            infile.write(line)
    os.remove('k_temp')
    os.remove('cl_temp')
    os.remove('mg_temp')
    os.remove('sol_temp')
    os.remove('biomol_temp')
    print('ok')
    print('\n')



# Writes the topology file
if make_top == 1:
# if [ $make_top -eq 1 ]; then
    print('\tWriting the topology file for the soup ... ')
    with open('temporary_top_part_1', 'w') as outfile:
        outfile.write('; Include forcefield parameters\n' )                  
        outfile.write('#include "../amber99sbws-dangK.ff/forcefield.itp"\n')
        outfile.write('#include "../gaff-types.itp"\n')                      
        outfile.write('\n')                                                
        outfile.write('; Include chain topologies\n')                       

	with open('temporary_top_part_3', 'w') as outfile:
        outfile.write('\n')
        outfile.write('; Include water topology\n')
        outfile.write('#include "../amber99sbws-dangK.ff/tip4p2005s.itp"\n')
        outfile.write('\n')
        outfile.write('; Include topology for ions\n' )  
        outfile.write('#include "../amber99sbws-dangK.ff/ions.itp"\n')
        outfile.write('\n')
        outfile.write('[ system ]\n')
        outfile.write('; Name\n')
        outfile.write('E.coli K12 cytoplasm model (aka soup)\n')
        outfile.write('\n')  
        outfile.write('[ molecules ]\n' )  
        outfile.write('; Compound        #mols\n' )
	
    with open('temporary_top_part_2', 'w') as outfile:
        for crowder in indices_sorted:
            outfile.write('#include "crowder_n%s.itp"' % (str(crowder)))
            outfile.write('crowder_n%s\t%s' % (str(crowder), str(num_crowder[crowder]) ))
    os.system('cp %s/crowder_n%s/crowder_n%s.itp %s' % (path_crowders, str(crowder),str(crowder), path_output))

    number_OW = float(sol_counter)/4
	with open('temporary_top_part_4', 'w') as outfile:
        outfile.write("SOL         %s" % str(number_OW))
        outfile.write("K         %s" % str(k_counter))
        outfile.write("CL         %s" % str(cl_counter))
        outfile.write("MG         %s" % str(mg_counter))
	
	os.system('cat temporary_top_part_1 temporary_top_part_2 temporary_top_part_3 temporary_top_part_4 > %s/top.top' % (path_output))
    os.remove('temporary_top_part_1')
    os.remove('temporary_top_part_2')
    os.remove('temporary_top_part_3')
    os.remove('temporary_top_part_4')
    print('ok') 
    print('\n')	



# if [ $make_neutral -eq 1 ]; then
if make_neutral == 1:
    os.system('%s grompp -f %s/%s -p %s/top.top -c %s/box_ordered.gro -o temporary.tpr .maxwarn 1 >out 2>err' % (args.gromacs_command, path_mdp, mdp_min, path_output))
    with open('err', 'r') as infile:
        for line in infile:
            if line.find('non-zero total charge') > -1:
                ions_neutralize = float(line.strip()[-1])
    os.remove('err')
    os.remove('out')
    os.remove('temporary.tpr')
    if ions_neutralize < 0:
        mass_ion_neutralize=39.1
        m1=-1
        ion_neutralize='K'
        ions_neutralize=ions_neutralize*m1
    else:
        mass_ion_neutralize=35
        ion_neutralize='Cl'
    print('\tAdding %s %s to neutralize the net charge of the system ... ' % (str(ions_neutralize), str(ion_neutralize)))
    os.system('%s grompp -f %s/%s -c %s/box_ordered.gro -p %s/top.top -o %s/ion_neutralize.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min, path_output, path_output, path_output))

    if ion_neutralize == 'K':
        os.system('echo SOL | %s genion -s %s/ion_neutralize.tpr -p %s/top.top -o %s/ion_neutralize.gro -np %s -pname %s >/dev/null 2>/dev/null' % (args.gromacs_command, path_output, path_output, ions_neutralize, ion_neutralize))
    else:
        os.system('echo SOL | %s genion -s %s/ion_neutralize.tpr -p %s/top.top -o %s/ion_neutralize.gro -nn %s -nname %s >/dev/null 2>/dev/null' % (args.gromacs_command, path_output, path_output, ions_neutralize, ion_neutralize))

    os.remove('mdout.mdp')
    print('ok')
    print('\n')


# Ionic Strength
if make_is == 1:
    with open('%s/ion_neutralize.gro' % (path_output), 'r') as infile:
        sol_counter=0
        for line in infile:
            if line.find('SOL'):
                sol_counter = sol_counter+1
    number_waters = sol_counter/4
    os.system('%s grompp -f %s/%s -c %s/ion_neutralize.gro -p %s/top.top -o %s/ion_is.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min, path_output, path_output, path_output))
    number_of_ions= ionic_strength*(number_waters/55.555)
    print('\tAdding %s KCl to achieve ionic strength of %s M ... ' % (str(number_of_ions), str(ionic_strength)))
    os.system('echo SOL | %s genion -s %s/ion_is.tpr -p %s/top.top -o %s/ion_is.gro -np %s -nn %s -pname K -nname CL >/dev/null 2>/dev/null' % (args.gromacs_command, path_output, path_output, number_of_ions, number_of_ions))
    os.remove('mdout.mdp')	
    print('ok')
    print('\n')

if check_fraction == 1:
    print('Check fraction')
# 	final_mass_water=$(grep ^"SOL " "$path_output"/top.top | awk '{print $2*18}')
# 	final_mass_MG=$(grep ^"MG " "$path_output"/top.top | awk '{sum +=$2} END {print sum*24.3}')
# 	final_mass_CL=$(grep ^"CL " "$path_output"/top.top | awk '{sum +=$2} END {print sum*35.4}')
# 	final_mass_K=$(grep ^"K " "$path_output"/top.top | awk '{sum +=$2} END {print sum*39.1}')
# 	total_nonwater_mass=$(echo "$final_mass_MG" "$final_mass_CL" "$final_mass_K" "$total_biomol_mass"| awk '{printf "%.2f\n", $1+$2+$3+$4}' )
	
# 	echo "final_mass_water: " "$final_mass_water"
# 	echo "final_mass_K: " "$final_mass_K"
# 	echo "final_mass_MG: " "$final_mass_MG"
# 	echo "final_mass_CL: " "$final_mass_CL"
# 	echo "total_nonwater_mass: " "$total_nonwater_mass"
	
# 	nonwater_mass_fraction=$(echo "$total_nonwater_mass" "$final_mass_water"  | awk '{printf "%.2f\n", $1/($1+$2)}' )
# 	biomol_fraction=$(echo $total_biomol_mass $final_mass_water $final_mass_MG $final_mass_CL $final_mass_K | awk '{printf "%.2f\n", $1/($1+$2+$3+$4+$5)}' )
	
# 	echo
# 	echo "nonwater     mass fraction = ""$nonwater_mass_fraction"
# 	echo "biomolecular mass fraction = ""$biomol_fraction"
# fi
