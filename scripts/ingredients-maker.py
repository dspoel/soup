#!/usr/bin/env python3

import argparse, os, glob, math
import subprocess
import numpy as np


path_mdp= "../mdp/" 


def parseArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gmx", "--gmx_path",  type=str,  default="gmx")
	parser.add_argument("-cr_n", "--n", help="crowder index",  type=int,  default=0)
	
	args = parser.parse_args()
	return args





def grp_list( gro ):
	# Creates a temporary index file to check what type of molecule is the crowder
#	print( 'q\n | gmx make_ndx -f crowder_n%d/npt.1.gro -o temporary.ndx >/dev/null 2>/dev/null ' %(crowder))
#	with open('q.txt', 'w') as outfile:
 #           outfile.write("q\n")
	os.system('echo "q""\\n" | %s make_ndx -f %s -o temporary.ndx >/dev/null 2>/dev/null' % (gmx_path, gro) )
	with open('temporary.ndx', 'r') as infile:
		lines = infile.readlines()
		grps = [x.split()[1] for x in lines if x[0]=='[']
	
	return grps
	os.remove('temporary.ndx')






def crowder_type( group_list ):
	if 'Protein' in group_list:
		return 'Protein', 'protein'
	elif 'RNA' in group_list:
		return 'RNA', 'nucleic acid'
	elif 'Other' in group_list:
		return 'Other', 'metabolite'
	else:
		return 'Unknown', 'unknown'






if __name__ == '__main__':
	args  = parseArguments()
	gmx_path = args.gmx_path
	crowder_num = args.n
	
	flag_build  = True
	flag_eq     = False
	flag_eq_nvt = False
	flag_eq_npt = False
	flag_prod   = False
	
	mdp_min_none = 'min_none.mdp'
	mdp_min      = 'min.mdp'
	mdp_nvt      = 'nvt.mdp'
	mdp_npt      = 'npt.mdp'
	mdp_prod     = 'prod.mdp'
	
	group_list = grp_list('crowder_n%d/crowder_n%d.pdb' % (crowder_num, crowder_num) )
	crowder_type , crowder_string = crowder_type( group_list )
	
	cation = {}
	cation['K'] = {}
	cation['K']['charge'] = 1
	cation['K']['conc'] = 0.150
	cation['MG'] = {}
	cation['MG']['charge'] = 2
	cation['MG']['conc'] = 0.06
	
	if crowder_string == 'protein':
		cation_type = 'K'
	elif crowder_string == 'nucleic acid':
		cation_type = 'MG'
	elif crowder_string == 'metabolite':
		if 'ATP' in group_list:
			cation_type = 'MG'
	else:
		cation_type = 'K'
	
	total_rep = 3
	
	
	
	if flag_build:
		print("\tCreating topolgy")	
		
		os.chdir('crowder_n%d' % crowder_num)
		
		with open('crowder_n%d.top' % crowder_num, 'w') as outfile:
			outfile.write('; Include forcefield parameters\n')
			outfile.write('#include "../amber99sbws-dangK.ff/forcefield.itp\n')
			if crowder_string == 'metabolite':
				outfile.write('#include "../gaff-types.itp\n')
			outfile.write('\n')
			outfile.write('; Include chain topologies\n')
			outfile.write('#include "crowder_n%s.itp\n' % crowder_num)
			outfile.write('\n')
			outfile.write('; Include water topology\n')
			outfile.write('#include "../amber99sbws-dangK.ff/tip4p2005s.itp\n')
			outfile.write('\n')
			outfile.write('; Include topologies for ions\n')
			outfile.write('#include "../amber99sbws-dangK.ff/ions.itp\n')
			outfile.write('\n')
			outfile.write('[ system ]\n')
			outfile.write('; Name\n')
			outfile.write('crowder_n%d\n' % crowder_num)
			outfile.write('\n')
			outfile.write('[ molecules ]\n')
			outfile.write('; Compound\t#mols\n')
			outfile.write('crowder_n%s\t1\n' % crowder_num)
		
		if crowder_string == 'metabolite':
			itps_list = glob.glob("*.itp")
			for itp in itps_list:
				words = itp.replace('.','_').split('_')
				if 'posre' in words:
					os.system('cp %s crowder_n%d_posres.itp' % (itp, crowder_num) )
					break
				else:
					os.system('cp %s temporary.top' % itp )
					break
		else:
			os.system('cp -r ../amber99sbws-dangK.ff/ .')
			
			os.system('echo 1 | %s pdb2gmx -f crowder_n%d.pdb -o crowder_n%d.gro -p temporary.top -i crowder_n%d_posres.itp -ignh -ff amber99sbws-dangK -merge all >/dev/null 2>/dev/null' % (args.gmx_path, crowder_num, crowder_num, crowder_num) )
			
			os.system('rm -r amber99sbws-dangK.ff/')
			
			
			
		output = False
		with open('temporary.top', 'r') as infile, open('crowder_n%d.itp' % crowder_num, 'w') as outfile:
			lines = infile.readlines()
			for line in lines:
				if line.strip() == '[ moleculetype ]':
					output = True
					n=1
				elif line.strip() == '; Include water topology':
					output = False
				
				if output:
					if n == 3:
						outfile.write('crowder_n%d         3\n' % crowder_num)
					else:
						outfile.write(line)
					n += 1
					
			
		os.remove('temporary.top')
		
		print('\tPutting the crowder in a box')
		os.system('echo System | %s editconf -f crowder_n%d.gro -o box.gro -c -d 2.0 -bt dodecahedron -princ >/dev/null 2>/dev/null' % (args.gmx_path, crowder_num))
		
		print('\tAdding water to the box')
		os.system('%s solvate -cp box.gro -cs ../amber99sbws-dangK.ff/tip4p2005.gro -o solv.gro -p crowder_n%d.top >/dev/null 2>/dev/null' % (args.gmx_path, crowder_num))
		
		print('\tAdding counter-ions to neutralize the charges')
		os.system('%s grompp -f %s/%s -c solv.gro -p crowder_n%d.top -o ions_neutralize.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gmx_path, path_mdp, mdp_min, crowder_num))
		
		os.system('echo SOL | %s genion -s ions_neutralize.tpr -p crowder_n%d.top -neutral -pname %s -pq %d -nname CL -nq -1 -o neutral.gro >/dev/null 2>/dev/null' % (args.gmx_path, crowder_num, cation_type, cation[cation_type]['charge']))
		
		print('\tAdding ions to model physiological ionic strenght')
		os.system('%s grompp -f %s/%s -c neutral.gro -p crowder_n%d.top -o ions_is.tpr >/dev/null 2>/dev/null' % (args.gmx_path, path_mdp, mdp_min, crowder_num))            
		os.system('echo SOL | %s genion -s ions_is.tpr -p crowder_n%d.top -conc %f -neutral -pname %s -pq %d -nname CL -nq -1 -o is.gro >/dev/null 2>/dev/null' % (args.gmx_path, crowder_num, cation[cation_type]['conc'], cation_type, cation[cation_type]['charge']) )
		
		print('\tEnergy minimization without constraints')
		os.system('%s grompp -f %s/%s -c is.gro -p crowder_n%d.top -o min_none.tpr >/dev/null 2>/dev/null' % (args.gmx_path, path_mdp, mdp_min_none, crowder_num))
		os.system('%s mdrun -v -pin on -deffnm min_none' % args.gmx_path)
		
		print('\tEnergy minimization with constraints')
		os.system('%s grompp -f %s/%s -c min_none.gro -p crowder_n%d.top -o min.tpr >/dev/null 2>/dev/null' % (args.gmx_path, path_mdp, mdp_min, crowder_num))
		os.system('%s mdrun -v -pin on -deffnm min' % args.gmx_path)
		
		os.chdir('..')





	if flag_eq:
		os.chdir('crowder_n%d' % crowder_num)
		
		# creates an index file with the group 'solute'. This is a way to call proteins, RNAs and metabolites by the same name from now on
		group_list = grp_list('min.gro')
		os.system('echo \'1\'"\n"\'name %d solute\'"\n"\'q\'"\n" | %s make_ndx -f min.gro -o index.ndx >/dev/null 2>/dev/null' % (len(group_list), args.gmx_path))
		
		
		for rep in range(1, total_rep+1):
			if flag_eq_nvt:
				if not os.path.isfile('nvt.%d.tpr' % rep):
					os.system('%s grompp -f %s/%s -c min.gro -r min.gro -n index.ndx -p crowder_n%d.top -o nvt.%d.tpr' % (args.gmx_path, path_mdp, mdp_nvt, crowder_num, rep))
				if not os.path.isfile('nvt.%d.gro' % rep):
					os.system('%s mdrun -s nvt.%d.tpr -v -pin on -deffnm nvt.%d -cpi nvt.%d.cpt -append' % (args.gmx_path, rep, rep, rep))
				
			if flag_eq_npt:
				if not os.path.isfile('npt.%d.tpr' % rep):
					os.system('%s grompp -f %s/%s -c nvt.%d.gro -r nvt.%d.gro -n index.ndx -p crowder_n%d.top -o npt.%d.tpr' % (args.gmx_path, path_mdp, mdp_npt, rep, rep, crowder_num, rep))
				if not os.path.isfile('npt.%d.gro' % rep):
					os.system('%s mdrun -s npt.%d.tpr -v -pin on -deffnm npt.%s -cpi npt.%d.cpt -append' % (args.gmx_path, rep, rep, rep))

		os.chdir('..')


	if flag_prod:
		os.chdir('crowder_n%d' % crowder_num)
		
		group_list = grp_list('npt.1.gro')
		os.system('echo \'1\'"\n"\'name %d solute\'"\n"\'q\'"\n" | %s make_ndx -f npt.1.gro -o index.ndx >/dev/null 2>/dev/null' % (len(group_list), args.gmx_path))
		
		for rep in range(1, total_rep+1):
			print('\tProduction run for replica number %d' % rep)
			if not os.path.isfile('prod.%d.tpr' % rep):
				os.system('%s grompp -f %s/%s -c npt.%d.gro -r nvt.%d.gro -n index.ndx -p crowder_n%d.top -o prod.%d.tpr' % (args.gmx_path, path_mdp, mdp_prod, rep, rep, crowder_num, rep))
			if not os.path.isfile('prod.%d.gro' % rep):
				os.system('%s mdrun -s prod.%s.tpr -v -pin on -deffnm prod.%s' % (args.gmx_path, rep, rep))
			
		os.chdir('..')
