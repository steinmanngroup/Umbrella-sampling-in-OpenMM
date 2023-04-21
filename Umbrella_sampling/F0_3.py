from tqdm import tqdm
from sys import argv
from datetime import datetime

import numpy as np

from FastMBAR.fastmbar import FastMBAR
from openmm.openmm import unit

from A0_setup import get_paths, get_par

start_time = datetime.now().timestamp()
script, main_path, run_on_gpu, high, L_times = argv
if run_on_gpu.lower() == 'yes':
    CUDA_run = True
    print('Utilizing the gpu where possible\n')
else:
    CUDA_run = False
    print('The FastMBAR will be a bit less fast, as it\'s running on the cpus instead..\n')

paths = get_paths(main_path)
distance_path = paths['distance_path']
graphs_path = paths['graphs_path']
data_path = paths['data_path']

# Force parameters
K_pull_prod = get_par('K_pull_prod', main_path)

# Integrator parameters
temp = get_par('temp', main_path)

# Umbrella sampling parameters
M = get_par('M', main_path)
R_min = get_par('R_min', main_path, get_only_value=True)
R_max = get_par('R_max', main_path, get_only_value=True)

Rs: list[float] = []
num_conf: list = []

if high.lower() in ['yes', 'y', 'ok']:
    high_interval = True
else:
    high_interval = False

if high_interval:
    M = int(M/2)
    R0:list = np.linspace(0, R_max, M, endpoint=False)
    L: int = M * int(L_times)
    data_file = f'{data_path}/data_hi_L{L}.csv'
    print('Reading the distances from the umbrella sampling')
    for R0_index in tqdm(range(M, M*2), colour='red'):
        R = np.loadtxt(f'{distance_path}/dist_{R0_index:03d}.csv')
        Rs.append(R)
        num_conf.append(len(R))
else:
    M = int(M/2)+1
    R0:list = np.linspace(abs(R_min), 0, M, endpoint=True)
    R0 = np.flip(R0)
    L: int = M * int(L_times)
    # data_file = f'{data_path}/data_lo_L{L}.csv'
    data_file = f'{data_path}/data_lo_reverse_L{L}.csv'
    print('Reading the distances from the umbrella sampling')
    for R0_index in tqdm(range(M), colour='red'):
        R = np.loadtxt(f'{distance_path}/dist_{R0_index:03d}.csv')
        Rs.append(R)
        num_conf.append(len(R))

print(f'\nThe interval is set to {R0[0]} to {R0[-1]} with {M} windows and {L} bins')
print('_'*100, '\n')
print(R0)

Rs = np.concatenate(Rs)
num_conf:list = np.array(num_conf).astype(np.float64)
N: int = len(Rs)
A: np.matrix = np.zeros([M, N])

kbT = unit.BOLTZMANN_CONSTANT_kB * temp * unit.AVOGADRO_CONSTANT_NA
kbt = kbT.value_in_unit(unit.kilocalorie_per_mole)
print('Creating the reduced energy matrix A')
for R0_index in tqdm(range(M), colour='green'):
    current_R0 = R0[R0_index]
    diff = np.abs(Rs-current_R0)
    print(current_R0, diff[0])
    A[R0_index, :] = 0.5 * K_pull_prod * diff ** 2/kbT



'''un_k is the reduced energy matrix (A), while N_k is the number of 
accepted conformations'''
start_time_fastmbar = datetime.now().timestamp()
fastmbar = FastMBAR(energy=A, num_conf=num_conf, cuda=CUDA_run, verbose=True, bootstrap=True, cuda_batch_mode=CUDA_run)

relative_G = fastmbar.F

end_time_fastmbar = datetime.now().timestamp()

print(f'\nRelative free energies: {relative_G}\n')

print('='*100)
print('Beginning the expansion of data now...')
print('-'*100)


R_PMF = np.linspace(abs(R0[0]), abs(R0[-1]), L, endpoint=False)
width = 10 / L
B = np.zeros((L, N))

for i in range(L):
    R_center = R_PMF[i]
    R_low = R_center - 0.5 * width
    R_high = R_center + 0.5 * width

    indicator = ((Rs > R_low) & (Rs <= R_high)) | \
        ((Rs + 1 > R_low) & (Rs + 1 <= R_high)) | \
            ((Rs - 1 > R_low) & (Rs - 1 <= R_high))

    B[i, ~indicator] = np.inf

PMF, std_error = fastmbar.calculate_free_energies_of_perturbed_states(B)
if CUDA_run:
    std_error = std_error.cpu().numpy()
print(f'\nSize of L vs M: {L}/{M} = {L/M}\n')
print('std error',std_error)
print('PMF', PMF)
min_PMF = np.min(PMF)
print('min_pmf:', min_PMF)


data = []
data.append(f'{"R0":<3s}\t{"PMF":<10s}\t{"ST_ERR":<7s}\t\n')
for i in range(len(PMF)):
    data.append(f'{R_PMF[i]:>3.2f}\t{PMF[i]:>10.2f}\t{std_error[i]:>7.2f}\n')


with open(data_file, 'w') as d_file:
    d_file.writelines(data)

if high == 'yes':
    print('Finished high now')
else:
    print('Finished low now')
end_time = datetime.now().timestamp()

print('#'*100)
print(f'\nIt took FastMBAR {datetime.fromtimestamp(abs(start_time_fastmbar-end_time_fastmbar))} to run')
print(f'\nThe whole program took {datetime.fromtimestamp(abs(start_time - end_time))}\n')
print('#'*100)