from sys import argv, stdout

from datetime import datetime

from openmm.app import PME
from openmm.app import Simulation
from openmm.app import PDBReporter
from openmm.app import StateDataReporter

from openmm.openmm import unit
from openmm.openmm import LangevinMiddleIntegrator
from openmm.openmm import MonteCarloBarostat

# Own functions
from create_CCBForce_system import get_SMIRNOFF_molecule_forcefield, create_modeller_box, add_pull_force_system
from alter_PDB import get_distance, get_host_serial, get_guest_serial
from setup import get_par, get_paths, get_ini_molecule_path

start_time = datetime.now().timestamp()
script, main_path, R0, R0_index = argv
R0_index = int(R0_index)
R0 = float(R0) * unit.angstroms

paths = get_paths(main_path)

initial_path = paths['initial_path']
minimize_path = paths['minimize_path']
displaced_path = paths['displace_path']
distance_path = paths['distance_path']
production_pdbs_path = paths['production_pdbs_path']
simulation_log_path = paths['simulation_log_path']



temp = get_par('temp', main_path)
fric_coef = get_par('fric_coef', main_path)
step_size = get_par('step_size', main_path)
pressure = get_par('pressure', main_path)

sdf_host = get_ini_molecule_path('host', main_path)
sdf_guest = get_ini_molecule_path('guest', main_path)
sdf_files = [sdf_host, sdf_guest]
forcefield = get_SMIRNOFF_molecule_forcefield(sdf_files)

print('Minimizing, equilibrating and producing the simulations.\n')

minimized_pdb = f'{minimize_path}/{R0_index:03d}.pdb'

host_serial = get_host_serial(minimized_pdb)
guest_serial = get_guest_serial(minimized_pdb)
    
# Setting up the moddeller and system with add_force_pull
modeller = create_modeller_box(in_pdb_file=minimized_pdb, forcefield=forcefield, water=True, box_side=5.0)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)

system = add_pull_force_system(host_serial, guest_serial, system)
print('Adding barostat')

barostat = MonteCarloBarostat(pressure, temp)
system.addForce(barostat)
print(f'The system is periodic: {system.usesPeriodicBoundaryConditions()}')
integrator = LangevinMiddleIntegrator(temp, fric_coef, step_size)

simulation = Simulation(modeller.topology, system, integrator)
print(f'The system has box vectors of: {simulation.topology.getPeriodicBoxVectors()}')
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(temp)

simulation.context.setParameter('R0', R0)

K_pull_min_eq = get_par('K_pull_min_eq', main_path)
simulation.context.setParameter('K_pull', K_pull_min_eq)

production_pdb = f'{production_pdbs_path}/R0_loop{R0_index:03d}.pdb'
simulation.reporters.append(StateDataReporter(stdout, reportInterval=500, step=True, 
                                              volume=True, temperature=True, potentialEnergy=True, 
                                              kineticEnergy=True, totalEnergy=True, separator='\t'))

window_number = f'{R0_index}'

print(f'{"Minimizing window":27}{window_number}:')
min_time = datetime.now().timestamp()
min_step = get_par('min_step', main_path)
simulation.minimizeEnergy(maxIterations=min_step)
actual_min_time = datetime.now().timestamp() - min_time
print(f'\nThe minization took {datetime.fromtimestamp(actual_min_time)}\n')

print(f'{"Equilibration of window":27}{window_number}:')
eq_time = datetime.now().timestamp()
eq_step = get_par('eq_step', main_path)
simulation.step(eq_step)
actual_eq_time = datetime.now().timestamp() - eq_time
print(f'\nThe equlibration took {datetime.fromtimestamp(actual_eq_time)}\n')

K_pull_prod = get_par('K_pull_prod', main_path)
simulation.context.setParameter('K_pull', K_pull_prod)


dist_file_path = f'{distance_path}/dist_{R0_index:03d}.csv'
open(dist_file_path, 'w').close()
filehandle = open(dist_file_path, 'a')

sim_time = datetime.now().timestamp()

sim_step = get_par('sim_step', main_path)
sim_loop = get_par('sim_loop', main_path)
report_interval = get_par('report_interval', main_path)
print(f'Production run of window {R0_index}')

simulation.reporters.append(StateDataReporter(stdout, reportInterval=report_interval, step=True, 
                                              volume=True, temperature=True, potentialEnergy=True, 
                                              kineticEnergy=True, totalEnergy=True, separator='\t'))
for i in range(sim_loop):
    production_pdb = f'{production_pdbs_path}/R0_loop{R0_index:03d}.pdb'
    simulation.reporters.append(PDBReporter(file=production_pdb, reportInterval=report_interval))
    
    
    simulation.step(sim_step)
    distance = get_distance(production_pdb)
    filehandle.write(f'{distance:.5f}\n')
    
    simulation.reporters.pop() # removes the last entry (PDBReporter but leaves StateDataReporters!!)

actuaL_sim_time = datetime.now().timestamp() - sim_time
print(f'\nThe simulation took {datetime.fromtimestamp(actuaL_sim_time)}\n')
actual_tot_time = datetime.now().timestamp() - start_time
print(f'\nThe whole program took {datetime.fromtimestamp(actual_tot_time)}\n')

# subprocess.run(['tar cvf'])

time_str = ['Minimization', 'Eqilibration', 'Production', 'Total']
time_log = [actual_min_time, actual_eq_time, actuaL_sim_time, actual_tot_time]

with open('time_report.txt', 'w') as t_report:
    header = f'{"Operation":15}\t{"time_stamp":10}\t{"time"}\n'
    t_report.write(header)
    print(header)
    basis_time = datetime.fromtimestamp(0.0)
    for i in range(len(time_str)):
        readable_time = str(abs(datetime.fromtimestamp(time_log[i])-basis_time))      
        line = f'{time_str[i]:15}\t{time_log[i]:10.3f}\t{readable_time}\n'
        print(line)
        t_report.write(line)


