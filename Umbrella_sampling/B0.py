from sys import argv

from openmm.app import NoCutoff
from openmm.openmm import LangevinMiddleIntegrator

# Own stuff
from create_CCBForce_system import get_SMIRNOFF_molecule_forcefield, create_modeller
from alter_PDB import write_centered_complex_to_orego, get_TER_pdb, get_centered_pdb
from setup import get_par, get_ini_molecule_path, get_paths
from create_CCBForce_system import write_system_to_file
script, main_path = argv

paths = get_paths(main_path)
main_path = paths['main_path']
initial_path = paths['initial_path']
minimize_path = paths['minimize_path']
displace_path = paths['displace_path']


sdf_host = get_ini_molecule_path('host', main_path)
sdf_guest = get_ini_molecule_path('guest', main_path)

sdf_files = [sdf_host, sdf_guest]
pdb_file_initial = get_ini_molecule_path('complex', main_path)

# Setting up the forceforcefield with SMIRNOFF
forcefield = get_SMIRNOFF_molecule_forcefield(sdf_files)
modeller = create_modeller(in_pdb_file=pdb_file_initial, forcefield=forcefield, water=False)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

write_system_to_file(system_variable=system, write_to_folder=initial_path)

# Integrator and variables

temp = get_par('temp', main_path)
fric_coef = get_par('fric_coef', main_path)
step_size = get_par('step_size', main_path)
print(f'The integrator is setup with:')
print(f'\t{"Temperature:":20}{temp} \n\t{"Step size:":20}{step_size} \n\t{"Friction coeffient:":20}{fric_coef}')

integrator = LangevinMiddleIntegrator(temp, fric_coef, step_size)
simulation_parameters = [modeller, system, integrator]
print('\nSetting the PDBFile up with terminating lines (TER)\n')
pdb_file_ter = get_TER_pdb(pdb_file_initial, initial_path, simulation_parameters)

print('Centering the host and guest to the same point in space\n')
integrator = LangevinMiddleIntegrator(temp, fric_coef, step_size)
simulation_parameters = [modeller, system, integrator]
pdb_file_center = get_centered_pdb(pdb_file_ter, initial_path, simulation_parameters)

print('\nMoving the complex to orego to make for a good starting point\n')
pdb_at_orego = write_centered_complex_to_orego(pdb_file_center, initial_path)