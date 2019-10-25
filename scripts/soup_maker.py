#!/usr/bin/env python3

import argparse, os, glob, math
import subprocess
import numpy as np

os.system('export GMX_MAXBACKUP=-1')


def parseArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gmx", "--gmx_path",  type=str,  default="gmx")
    
	args = parser.parse_args()
	return args




box_size_initial=35 #35		# in nm
box_size_step=2	# nm

droplet_thickness_initial=0.6 #0.6
droplet_shell_precision=2	# percentage as integer

crowder_total=13
biomolecular_fraction=0.30
ionic_strength=0.150	# M

make_droplets = False
make_box      = False
make_topology = True

path_mdp='mdp'
path_output='soup'

mdp_min='min.mdp'



# number of gro files for each crowder


crowders = {}
for c in range(1,crowder_total+1):
	crowders[c] = {}

crowders[1]['copies']  = 6
crowders[2]['copies']  = 7
crowders[3]['copies']  = 2
crowders[4]['copies']  = 1
crowders[5]['copies']  = 3
crowders[6]['copies']  = 1
crowders[7]['copies']  = 1
crowders[8]['copies']  = 1
crowders[9]['copies']  = 5
crowders[10]['copies'] = 1436
crowders[11]['copies'] = 144
crowders[12]['copies'] = 225
crowders[13]['copies'] = 255




def grp_list( gro ):
	# Creates a temporary index file to check what type of molecule is the crowder
#	print( 'q\n | gmx make_ndx -f crowder_n%s/npt.1.gro -o temporary.ndx >/dev/null 2>/dev/null ' %(crowder))
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

def crowder_mass( crowder_id ):
	with open('crowder_n%d/crowder_n%d.itp' % (crowder_id, crowder_id), 'r') as infile:
		mass = np.array([float(x.split()[-4]) for x in infile.readlines() if (('qtot' in x.split()) and (x.split()[0]!=';'))]).sum()
	return mass

def count_waters( gro_file ):
	waters = 0
	with open(gro_file, 'r') as infile:
		for line in infile.readlines():
			try:
				words = line.split()
				if (words[0][-3:]=='SOL') and (words[1][:2]=='OW'):
					waters += 1
			except:
				pass
	return waters
            





if __name__=='__main__':
	
	args  = parseArguments()
	gmx_path = args.gmx_path
	
	
	print('Gathering data about each crowder')
	
	
	for crowder_id in range(1, crowder_total+1):
		# Creates a temporary index file to check what type of molecule is the crowder
		grps = grp_list('crowder_n'+str(crowder_id)+'/npt.1.gro')
		crowders[crowder_id]['type'], crowders[crowder_id]['class'] = crowder_type(grps)
		crowders[crowder_id]['mass'] = crowder_mass(crowder_id)
		
		print('\tCrowder %d is a %s with mass = %.2f Da' % (crowder_id, crowders[crowder_id]['class'], crowders[crowder_id]['mass']) )
	
	
	
	biomolecular_mass=0
	for crowder_id in crowders:
		biomolecular_mass += crowders[crowder_id]['copies'] * crowders[crowder_id]['mass']

	print('')
	print('\tTotal biomolecular mass = %.2f Da' % (biomolecular_mass))
	print('\tCalculating the number of waters necessary to achieve the biomolecular fraction of %f %%' % (biomolecular_fraction*100))
	
	water_mass = (biomolecular_mass/biomolecular_fraction) - biomolecular_mass
	water_molecules = int(water_mass / 18)
	
	print('\t%d water molecules are necessary to reach the biomolecular fraction of %.1f %%' % (water_molecules, biomolecular_fraction*100))
	print('')
	print('')
	
	crowders_sorted = sorted( crowders.items() , key= lambda tup: tup[1]['mass'], reverse=True)
	crowders_order = [x[0] for x in crowders_sorted]
	
	if make_droplets:
		print('Making droplets for each crowder')
		print('\tFinding the optimum droplet shell thickness which comprises a total of %d water molecules for all crowders\n' % water_molecules)
		print('\tThe maximum difference that is allowed between the number of water molecules inside the droplets and the number of molecules we need is %d %%' % droplet_shell_precision)
		
		droplet_thickness_step=0.1	# nm
		droplet_thickness=droplet_thickness_initial		# nm
		droplet_thickness_prev=droplet_thickness
		
		
		while True:
			water_added_to_the_soup = 0
			
			for crowder_id in range(1,crowder_total+1):
				print('\t\tMaking droplet for crowder %d ... ' % crowder_id, end='')
				
				os.chdir('crowder_n%d' % crowder_id)
				
				os.system('echo "System" | %s trjconv -f npt.1.gro -s npt.1.tpr -pbc atom -ur compact -o eq_pbc.gro >/dev/null 2>/dev/null' % args.gmx_path)
				
				grps = grp_list('npt.1.gro')
				
				if 'MG' in grps:
					grp_mg = 0
					grp_count = 0
					for grp in grps:
						if (grp == 'MG') and (grp_mg == 0):
							grp_mg = grp_count
						grp_count += 1
						
					os.system('%s select -f eq_pbc.gro -on index_droplet.ndx -select "(group %s ) or (same residue as (group "SOL" and within %f of group %s) ) or (group %d and within %s of group %s )" -s npt.1.tpr >/dev/null 2>/dev/null' % (args.gmx_path, crowders[crowder_id]['type'], droplet_thickness, crowders[crowder_id]['type'], grp_mg, droplet_thickness, crowders[crowder_id]['type']))
				else:
					os.system('%s select -f eq_pbc.gro -on index_droplet.ndx -select "(group %s ) or (same residue as (group "SOL" and within %f of group %s) )" -s npt.1.tpr >/dev/null 2>/dev/null' % (args.gmx_path, crowders[crowder_id]['type'], droplet_thickness, crowders[crowder_id]['type']))
				
				os.system('%s trjconv -f eq_pbc.gro -s npt.1.tpr -o droplet.gro -n index_droplet.ndx >/dev/null 2>/dev/null' % args.gmx_path)
				
				
				water_added_by_one_droplet = count_waters('droplet.gro')
				
				water_added_by_this_crowder = water_added_by_one_droplet * crowders[crowder_id]['copies']
				
				water_added_to_the_soup += water_added_by_this_crowder
				
				os.chdir('..')
				print('ok')
			
			water_added_perc = float(water_added_to_the_soup/water_molecules)*100
			error_perc = np.absolute(water_added_perc - 100)
			
			current_biomolecular_fraction = biomolecular_mass / (biomolecular_mass + (water_added_to_the_soup*18.0))
			
			print('\t\t\tDroplet thickness = %f' % droplet_thickness)
			print('\t\t\tCurrent biomolecular fraction = %f' % current_biomolecular_fraction)
			print('\t\t\t%% of water added = %f' % water_added_perc)
			print('\t\t\t%% error = %f' % error_perc)
			
			if error_perc <= droplet_shell_precision:
				break
			else:
				if water_added_perc < 100 :
					droplet_thickness_prev = droplet_thickness
				else:
					droplet_thickness = droplet_thickness_prev
					droplet_thickness_step = droplet_thickness_step/10
				
				droplet_thickness = droplet_thickness + droplet_thickness_step
		
		print('')
		print('\tDroplet shell thickness converged to %.1f nm.' % droplet_thickness)
		print('\t%d water molecules will be added to the box.' % water_added_to_the_soup)
		print('\tThis is %.2f %% of the precise number required' % water_added_perc)
		print('\tThat\'s within the requested error margin of %.2f %%' % droplet_shell_precision)
		print('')
		print('')



	if make_box:
		print('Finding the smallest cubic box in which all components of the soup fit')
		
		print('\tSorting the crowders by mass')
		
		for crowder_id in crowders_order:
			print('\t\tCrowder %d, mass = %d' % (crowder_id, crowders[crowder_id]['mass']) )
		print('')
		
		try:
			os.mkdir(path_output)
		except:
			pass
		
		box_size = box_size_initial
		ok = False
		while True:
			with open('%s/box.gro' % (path_output), 'w') as outfile:
				outfile.write('Cytoplasm model\n')
				outfile.write('0\n')
				outfile.write('%.4f\t%.4f\t%.4f\n' % (box_size, box_size, box_size))
			print('\tSide lenght = %.4f nm' % box_size)
			
			
			
			crowder_count = 1
			for crowder_id in crowders_order:
				print('\tTrying to add %d copies of crowder number %d ... ' % (crowders[crowder_id]['copies'], crowder_id), end='')
				os.system('%s insert-molecules -f %s/box.gro -ci crowder_n%d/droplet.gro -o %s/box.gro -nmol %d -try 10 >out 2>err' % (gmx_path, path_output, crowder_id, path_output, crowders[crowder_id]['copies']))
				with open('err', 'r') as infile:
					for line in infile.readlines():
						try:
							if line.split()[0]=='Added':
								number_added = int(line.split()[1])
						except:
							pass
				os.remove("err")
				os.remove("out")
				print('added %d ...' % number_added, end='')
				
				if number_added == crowders[crowder_id]['copies']:
					print('ok')
					if crowder_count == crowder_total:
						ok = True
						break
					else:
						crowder_count += 1
				else:
					print('didn\'t fit')
					box_size = box_size + box_size_step
					break
			if ok:
				break



if make_topology:
	print('')
	print('Writing the .gro file for the soup ... ', end='')
	
	sections = ['title', 'atoms', 'biomol', 'sol', 'mg', 'box']
	gro = {}
	for grp in sections:
		gro[grp] = []
	
	with open('%s/box.gro' % path_output , 'r') as infile:
		lines = infile.readlines()
		
		gro['title'].append(lines[0])
		gro['atoms'].append(lines[1])
		gro['box'].append(lines[-1])
		
		del lines[-1]
		del lines[:2]
		
		for line in lines:
			words = line.split()
			if words[0][-3:] == 'SOL':
				gro['sol'].append(line)
			elif words[0][-2:] == 'MG' and words[1][:2]=='MG':
				gro['mg'].append(line)
			else:
				gro['biomol'].append(line)
	
		with open('%s/box_ordered.gro' % path_output, 'w') as outfile:
			for grp in sections:
				for line in gro[grp]:
					outfile.write(line)
	print('ok')
		
		
	print('Writing the topology file for the soup ... ', end='')
	num_sol = int(len(gro['sol'])/4)
	num_mg  = int(len(gro['mg']))
	with open('%s/topology.top' % path_output, 'w') as outfile:
		outfile.write('; Include forcefield parameters\n' )                  
		outfile.write('#include "../amber99sbws-dangK.ff/forcefield.itp"\n')
		outfile.write('#include "../gaff-types.itp"\n')                      
		outfile.write('\n')                                                
		outfile.write('; Include chain topologies\n')
		for crowder_id in crowders_order:
			outfile.write('#include "crowder_n%d.itp"\n' % crowder_id)
			os.system('cp crowder_n%d/*.itp %s/' % (crowder_id, path_output) )
		outfile.write('\n')
		outfile.write('; Include water topology\n')
		outfile.write('#include "../amber99sbws-dangK.ff/tip4p2005s.itp"\n')
		outfile.write('\n')
		outfile.write('; Include topology for ions\n' )  
		outfile.write('#include "../amber99sbws-dangK.ff/ions.itp"\n')
		outfile.write('\n')
		outfile.write('[ system ]\n')
		outfile.write('; Name\n')
		outfile.write('E.coli K12 cytoplasm model\n')
		outfile.write('\n')  
		outfile.write('[ molecules ]\n' )  
		outfile.write('; Compound        #mols\n' )
		for crowder_id in crowders_order:
			outfile.write('crowder_n%d\t%d\n' % (crowder_id, crowders[crowder_id]['copies']) )
		outfile.write('SOL\t%d\n' % num_sol )
		outfile.write('MG\t%d\n' % num_mg )
	print('ok')
	print('')
	print('')
	
	
	
	print('Neutralizing charges in the box with K+ or Cl- ... ', end='')
	os.chdir(path_output)
	os.system('%s grompp -f ../mdp/%s -p topology.top -c box_ordered.gro -o ion_neutralize.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gmx_path, mdp_min))
	os.system('echo SOL | %s genion -s ion_neutralize.tpr -p topology.top -o ion_neutralize.gro -pname K -nname CL -pq 1 -nq -1 -neutral >/dev/null 2>/dev/null' % args.gmx_path)
	os.remove('mdout.mdp')
	print('ok')
	
	
	print('Adding KCl to reach the ionic strength of %f mol/L ... ' % ionic_strength , end='' )
	num_kcl = int(ionic_strength * (num_sol/55.555))
	os.system('%s grompp -f ../mdp/%s -p topology.top -c ion_neutralize.gro -o ion_ionicstrength.tpr -maxwarn 1 >/dev/null 2>/dev/null' % (args.gmx_path, mdp_min))
	os.system('echo SOL | %s genion -s ion_ionicstrength.tpr -p topology.top -o ion_ionicstrength.gro -pname K -nname CL -pq 1 -nq -1 -np %d -nn %d >/dev/null 2>/dev/null' % (args.gmx_path, num_kcl, num_kcl))
	os.remove('mdout.mdp')
	print('ok')
