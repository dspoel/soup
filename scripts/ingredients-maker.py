#!/usr/bin/env python3

import argparse, os, glob, math
import subprocess
import numpy as np

path="~/soup/"
path_mdp= "../mdp/" 
# path_mdp= "%s/mdp/" % (path)

def parseArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gmx", "--gromacs_command",  type=str,  default="gmx")
	parser.add_argument("-cr_n", "--n", help="crowder index",  type=int,  default=0)
	
	args = parser.parse_args()
	return args



# export GMX_MAXBACKUP=-1

if __name__ == '__main__':
    args  = parseArguments()
    gromacs_command = args.gromacs_command
    n=args.n 				#crowder index
    crowder_total=n
    rep_first=1
    rep_total=3
    flag_build=1
    flag_eq=0
    flag_eq_nvt=0
    flag_eq_npt=0
    flag_prod=1
    mdp_min_none='min_none.mdp'
    mdp_min='min.mdp'
    mdp_nvt='nvt.mdp'
    mdp_npt='npt.mdp'
    mdp_prod='prod.mdp'
    topology_template='%s/topology_template.top' % (path)
    while n <= crowder_total:
        print('Crowder No. %s:\n' % (str(n)))
        print('\tChecking the ingredient ... ')
        # Creates a temporary index file to check what type of molecule is the crowder
        with open('q.txt', 'w') as outfile:
            outfile.write("q\n")
        os.system('%s make_ndx -f crowder_n%s/%s_*pdb -o temporary.ndx < q.txt  >/dev/null 2>/dev/null ' % (gromacs_command, str(n), str(n) ))
        os.system('rm q.txt')
        with open('temporary.ndx', 'r') as infile:
            l=[]
            for line in infile:
                if line.find('[') > -1:
                    l.append(line.strip())
            if '[ Protein ]' in l:
                crowder_type = "Protein"
            elif '[ RNA ]' in l:
                crowder_type = "RNA"
            elif '[ Other ]' in l and '[ ATP ]' in l:
                crowder_type = "Other"
                is_atp = 1
            elif '[ Other ]' in l and not '[ ATP ]' in l:
                crowder_type = "Other"
                is_atp = 0
            else:
                print("Unknown type for crowder_n%s" % str(n))
        
        # At this point we have the crowder type stored in the crowder_type variable
        # Is is either 'Protein', 'RNA' or 'Other'
        # We also know if the metabolite is ATP

        print('it is a', end='')

        if crowder_type == "Protein":
          print(' protein')
          cation_type='K'
          cation_charge=1
          cation_conc=0.150
            
        elif crowder_type == "RNA":
          print('n RNA')
          cation_type='MG'
          cation_charge=2
          cation_conc=0.06
                
        elif crowder_type == "Other":
          print(' metabolite', end = '')
          if is_atp == 1:
              print(', and it is ATP')
              cation_type='MG'
              cation_charge=2
              cation_conc=0.06
          else:
              print(', and it is not ATP')
              cation_type='K'
              cation_charge=1
              cation_conc=0.150

        os.system('rm temporary.ndx')

        if flag_build == 1:
            print("\tCreating topolgy")
        #   # If the crowder is a metabolite, find the topology from gaff which is already available
            if crowder_type == "Other":
                # temporary_topology=$(find 'crowder_n'"$n"'/' -maxdepth 1 -name "*.itp" | grep -v 'posre')
                if glob.glob("crowder_n%s/*posre*itp" % (str(n))):
                    temporary_topology = glob.glob("crowder_n%s/*posre*itp" % (str(n)))[0]
                elif glob.glob("crowder_n%s/*posre*top" % (str(n))):
                    temporary_topology = glob.glob("crowder_n%s/*posre*top" % (str(n)))[0]
                else:
                    print("No topology file for crowder_n%s" % (str(n)))
                temporary_posre = glob.glob("crowder_n%s/*posre*itp" % (str(n)))[0]
              # in order to ensure naming consistency among all metabolites ...                    
                os.system("%s editconf -f crowder_n%s/%s_*pdb -o crowder_n%s/crowder_n%s.gro " % (args.gromacs_command, str(n), str(n), str(n), str(n)))
                os.system('cp temporary_topology crowder_n%s/temporary.top' % (str(n)))
                os.system('cp temporary_posre crowder_n%s/crowder_n%s_posre.itp 1>/dev/null' % (str(n), str(n)))
            else:
              # If the crowder is not a metabolite, we will use pdb2gmx to generate the top, itp and gro files
                os.system('echo 1 | %s pdb2gmx -f crowder_n%s/%s_*pdb -o crowder_n%s/crowder_n%s.gro -p crowder_n%s/temporary.top -i crowder_n%s/crowder_n%s_posres.itp -ignh -ff amber99sbws-dangK -merge all >/dev/null 2>/dev/null ' % (args.gromacs_command, str(n), str(n), str(n), str(n), str(n), str(n), str(n)))

            # extract relevant part of the itp/top file
            # this is done to ensure naming consistency
            with open('crowder_n%s/temporary.top' % (str(n)), 'r') as infile:
                text_list = infile.readlines()
                line_start = [i for i, elem in enumerate(text_list) if '[ atoms ]' in elem]
                line_end = [i for i, elem in enumerate(text_list) if '#ifdef POSRES\n' in elem]

            with open('temporary_body', 'w') as outfile:
                for item in text_list[line_start[0] : line_end[0]-2]:
                    outfile.write("%s\n" % item)

            with open('temporary_header', 'w') as outfile:
                outfile.write('[ moleculetype ]\n')
                outfile.write('; Name\t\tnrexcl\n')
                outfile.write('crowder_n%s\t\t3\n' % (str(n)))
                outfile.write('\n')
                
            with open('temporary_footer', 'w') as outfile:
                outfile.write('#ifdef POSRES\n')
                outfile.write('#include "crowder_n%s_posres.itp"\n' % (str(n)))
                outfile.write('#endif\n')

            os.system('cat temporary_header temporary_body temporary_footer > crowder_n%s/crowder_n%s.itp' % (str(n), str(n)))
            if os.path.isfile('crowder_n%s/crowder_n%s.itp' % (str(n), str(n))):
                os.system('rm temporary_header temporary_body temporary_footer')

            with open('crowder_n%s/crowder_n%s.top' % (str(n), str(n)), 'w') as outfile:
                outfile.write('; Include forcefield parameters\n')
                outfile.write('#include "../amber99sbws-dangK.ff/forcefield.itp\n')
                outfile.write('#include "../gaff-types.itp\n')
                outfile.write('\n')
                outfile.write('; Include chain topologies\n')
                outfile.write('#include "crowder_n%s.itp\n' % (str(n)))
                outfile.write('\n')
                outfile.write('; Include water topology\n')
                outfile.write('#include "../amber99sbws-dangK.ff/tip4p2005s.itp\n')
                outfile.write('\n')
                outfile.write('#ifdef POSRES_WATER\n')
                outfile.write('; Positin restraint for each water oxygen\n')
                outfile.write('[ position_restraints ]\n')
                outfile.write('; i\tfunct\tfcs\tfcy\tfcz\n')
                outfile.write('1\t1\t1000\t1000\t1000\n')
                outfile.write('#endif\n')
                outfile.write('\n')

                outfile.write('; Include topologies for ions\n')
                outfile.write('#include "../amber99sbws-dangK.ff/ions.itp\n')
                outfile.write('\n')
                outfile.write('[ system ]\n')
                outfile.write('; Name\n')
                outfile.write('crowder_n%s\n' % (str(n)))
                outfile.write('\n')
                outfile.write('[ molecules ]\n')
                outfile.write('; Compound\t#mols\n')
                outfile.write('crowder_n%s\t1\n' % (str(n)))


            os.chdir('crowder_n%s/' % (str(n)))
            print('\tPutting the crowder in a box')
            os.system('echo System | %s editconf -f crowder_n%s.gro -o box.gro -c -d 2.2 -bt dodecahedron -princ >/dev/null 2>/dev/null' % (args.gromacs_command, str(n)))

            print('\tAdding water to the box')
            os.system('%s solvate -cp box.gro -cs ../amber99sbws-dangK.ff/tip4p2005.gro -o solv.gro -p crowder_n%s.top >/dev/null 2>/dev/null' % (args.gromacs_command, str(n)))
            os.system('%s grompp -f %s/%s -c solv.gro -p crowder_n%s.top -o select.tpr -maxwarn 1 ' % (args.gromacs_command, path_mdp, mdp_min, str(n)))

            print('\tRemoving a water shell around the protein')
            os.system('%s select -f solv.gro -on index_solv.ndx -select "Close to protein" not (same residue as((resname SOL) and (within 0.5 of group "non-Water"))) -s select.tpr >/dev/null 2>/dev/null' % (args.gromacs_command))
            os.system('%s trjconv -f solv.gro -s select.tpr -o removed_water_shell.gro -n index_solv.ndx' % (args.gromacs_command))
            os.system('rm index_solv.ndx select.tpr')

            with open('removed_water_shell.gro', 'r') as infile:
                text = infile.readline()
                water_number = text.count("OW")
            with open('crowder_n%s.top' % (str(n)), 'r') as infile:
                with open('temp.top', 'w') as outfile:
                    for line in infile:
                        if line.find('SOL') < 0:
                            outfile.write(line)
            os.system('cp temp.top crowder_n%s.top' % (str(n)))
            with open('crowder_n%s.top' % (str(n)), 'a') as outfile:
                outfile.write('SOL\t%s\n' % (str(water_number)))
            
            print('\tAdding counter-ions to neutralize the charges')
            os.system('%s grompp -f %s/%s -c removed_water_shell.gro -p crowder_n%s.top -o ions_neutralize.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min, str(n)))
            os.system('echo \t SOL | %s genion -s ions_is.tpr -p crowder_n%s.top -conc 0.000 -neutral -pname %s -pq %s -nname CL -o neutral.gro >/dev/null 2>/dev/null' % (args.gromacs_command, str(n), str(cation_type), str(cation_charge)))
            
            print('\tAdding ions to model physiological ionic strenght')
            os.system('%s grompp -f %s/%s -c neutral.gro -p crowder_n%s.top -o ions_is.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min, str(n)))            
            os.system('echo \t SOL | %s genion -s ions_is.tpr -p crowder_n%s.top -conc %s -neutral -pname %s -pq %s -nname CL -o is.gro >/dev/null 2>/dev/null' % (args.gromacs_command, str(n), str(cation_conc), str(cation_type), str(cation_charge)))

            print('\tEnergy minimization without constraints')
            os.system('%s grompp -f %s/%s -c is.gro -p crowder_n%s.top -p min_none.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min_none, str(n)))
            os.system('%s mdrun -v -pin on -deffnm min_non >/dev/null 2>/dev/null')

            print('\tEnergy minimization with constraints')
            os.system('%s grompp -f %s/%s -c min_none.gro -p crowder_n%s.top -p min.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_min, str(n)))
            os.system('%s mdrun -v -pin on -deffnm min >/dev/null 2>/dev/null')

            os.chdir('..')

        if flag_eq == 1:
            os.chdir('crowder_n%s' % (str(n)))
          # creates an index file with the group 'solute'. This is a way to call proteins, RNAs and metabolites by the same name from now on
            with open('q.txt', 'w') as outfile:
                otufile.write("q\n")
            os.system('%s make_ndx -f min.gro -o index.ndx < q.txt >/dev/null 2>/dev/null')
            with open('index.ndx', 'r') as infile:
                file = infile.readlines()
                number_group = file.count('[')
            with open('q.txt', 'w') as outfile:
                outfile.write('1\n')
                outfile.write('name\t%s solute\n' % (str(number_group)))
                otufile.write("q\n")
            os.system('%s make_ndx -f min.gro -o index.ndx < q.txt >/dev/null 2>/dev/null ')

            rep=rep_first
            while rep <= rep_total:
                print('\tEquilibration/Thermalization run for replica number %s' % (str(rep)))
                # NVT equilibration
                if flag_eq_nvt == 1:
                    if not os.path.isfile('nvt.%s.tpr' % (str(rep))):
                        os.system('%s grompp -f %s/%s -c min.gro -r min.gro -n index.ndx -p crowder_n%s.top -o nvt.%s.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_nvt, str(n), str(rep)))
                    if not os.path.isfile('nvt.%s.gro' % (str(rep))):
                        os.system('%s mdrun -s nvt.%s.tpr -v -pin on -deffnm nvt.%s -cpi nvt.%s.cpt -append  >/dev/null 2>/dev/null' % (args.gromacs_command, str(rep), str(rep), str(rep)))
                # NPT equilibration
                if flag_eq_npt == 1:
                    if not os.path.isfile('npt.%s.tpr' % (str(rep))):
                        os.system('%s grompp -f %s/%s -c nvt.%s.gro -r nvt.%s.gro -n index.ndx -p crowder_n%s.top -o npt.%s.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_npt, str(rep), str(rep), str(n), str(rep)))
                    if not os.path.isfile('npt.%s.gro' % (str(rep))):
                        os.system('%s mdrun -s npt.%s.tpr -v -pin on -deffnm npt.%s -cpi npt.%s.cpt -append  >/dev/null 2>/dev/null' % (args.gromacs_command, str(rep), str(rep), str(rep)))
                rep = rep+1
            os.chdir('..')

        if flag_prod == 1:
            os.chdir('crowder_n%s' % (str(n)))
            # creates an index file with the group 'solute'. This is a way to call proteins, RNAs and metabolites by the same name from now on
            with open('q.txt', 'w') as outfile:
                otufile.write("q\n")
            os.system('%s make_ndx -f min.gro -o index.ndx < q.txt >/dev/null 2>/dev/null')
            with open('index.ndx', 'r') as infile:
                file = infile.readlines()
                number_group = file.count('[')
            with open('q.txt', 'w') as outfile:
                outfile.write('1\n')
                outfile.write('name\t%s solute\n' % (str(number_group)))
                otufile.write("q\n")
            os.system('%s make_ndx -f min.gro -o index.ndx < q.txt >/dev/null 2>/dev/null ')
            os.system('rm q.txt')
            while rep <= rep_total:
                print('\tProduction run for replica number %s' % (str(rep)))
                if not os.path.isfile('prod.%s.tpr' % (str(rep))):
                    os.system('%s grompp -f %s/%s -c npt.%s.gro -r nvt.%s.gro -n index.ndx -p crowder_n%s.top -o prod.%s.tpr >/dev/null 2>/dev/null' % (args.gromacs_command, path_mdp, mdp_prod, str(rep), str(rep), str(n), str(rep)))
                if not os.path.isfile('npt.%s.gro' % (str(rep))):
                    os.system('%s mdrun -s prod.%s.tpr -v -pin on -deffnm prod.%s   >/dev/null 2>/dev/null' % (args.gromacs_command, str(rep), str(rep)))
                rep = rep+1
        n=n+1











