from sys import argv

from openmm.app import Simulation
from openmm.app import PDBReporter
from openmm.openmm import LangevinMiddleIntegrator
from openmm.openmm import unit
 
# Own stuff
from create_CCBForce_system import create_modeller, get_SMIRNOFF_molecule_forcefield
from create_CCBForce_system import add_pull_force_system
from alter_PDB import displace_guest, get_distance, get_host_serial, get_guest_serial
from setup import get_par, get_paths, get_ini_molecule_path

script, main_path, R0, R0_index = argv
R0_index = int(R0_index)
R0 = float(R0)

paths = get_paths(main_path)

initial_path = paths['initial_path']
displace_path = paths['displace_path']
minimize_path = paths['minimize_path']
pdb_at_orego = f'{initial_path}/orego.pdb'


temp = get_par('temp', main_path)
fric_coef = get_par('fric_coef', main_path)
step_size = get_par('step_size', main_path)

print(f'Displacing the {R0_index}th complex with a factor of {R0}:')
displaced_pdb = displace_guest(pdb_at_orego, R0, R0_index, displace_path)
print(f'Distance after displacement of pdb: {get_distance(displaced_pdb):8.4f} angstrom')
print(displaced_pdb)

sdf_host = get_ini_molecule_path('host', main_path)
sdf_guest = get_ini_molecule_path('guest', main_path)
sdf_files = [sdf_host, sdf_guest]

forcefield = get_SMIRNOFF_molecule_forcefield(sdf_files)

modeller = create_modeller(displaced_pdb, forcefield=forcefield, water=False)
system = forcefield.createSystem(modeller.topology)

host_serial = get_host_serial(displaced_pdb)
guest_serial = get_guest_serial(displaced_pdb)
system = add_pull_force_system(host_serial, guest_serial, system)

integrator = LangevinMiddleIntegrator(temp, fric_coef, step_size)
simulation = Simulation(modeller.topology, system, integrator)


K_pull_min_eq = get_par('K_pull_min_eq', main_path)
simulation.context.setParameter('K_pull', K_pull_min_eq)

print('Minimizing the displaced structures...')

simulation.context.setPositions(modeller.positions)
simulation.context.setParameter('R0', abs(R0) * unit.angstrom)
displace_min_step = get_par('displace_min_step', main_path)
simulation.minimizeEnergy(maxIterations=displace_min_step)

minimized_pdb = f'{minimize_path}/{R0_index:03d}.pdb'
report_steps = 1
simulation.reporters.append(PDBReporter(minimized_pdb, reportInterval=report_steps))
simulation.step(report_steps)
print(f'dist = {get_distance(minimized_pdb):8.4f} after being minimized.')

print(
'''
json til dict - "det' nemt" - Casper Steinmann 30/11-2021 kl. 11.40
'''
)

print(f'The {R0_index}th complex has been displaced and minimized')