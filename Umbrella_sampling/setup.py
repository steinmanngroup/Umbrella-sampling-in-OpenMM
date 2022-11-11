import csv
from shutil import copy
from pathlib import Path
import os
import numpy as np
from collections import OrderedDict
from openmm.openmm import unit


def get_paths(main_path:str):
    '''
    Takes the main_path as a string and returns a dictionary.
        
    Setup:
        paths = get_paths('main_path')
        
        initial_path = paths['initial_path']
        
    Possible entries are:
    - 'main_path'
    - 'initial_path'
    - 'minimize_path'
    - 'displace_path'
    - 'production_path'
    - 'production_pdbs_path'
    - 'distance_path'
    - 'simulation_log_path'
    - 'graphs_path'
    - 'data_path'
    '''
    # Paths
    paths:dict[str, str] = {
        'main_path'             :   main_path,
        'initial_path'          :   f'{main_path}/start',
        'minimize_path'         :   f'{main_path}/minimize',
        'displace_path'         :   f'{main_path}/displace',
        'production_path'       :   f'{main_path}/production',
        'production_pdbs_path'  :   f'{main_path}/production/pdbs',
        'distance_path'         :   f'{main_path}/production/distance',
        'simulation_log_path'   :   f'{main_path}/production/logs',
        'graphs_path'           :   f'{main_path}/graphs',
        'data_path'             :   f'{main_path}/data'
    }
    return paths

def set_paths(main_path:str):
    '''
    Setup the pathing for the simulations.
        Provide the mainpath as a string.
        
    In case the folders already exists nothing is deleted.
    '''
    # cover_folder = f'{main_path}/{main_path}'
    paths = get_paths(main_path)
    
    for path_key in paths.keys():
        make_path = f'{paths[f"{path_key}"]}'
        Path(make_path).mkdir(parents=True, exist_ok=True)





def read_parameters(setup_file:str):
    with open(setup_file, 'r') as s_file:
        par_reader = csv.reader(s_file, delimiter='\t')
        par_rows =[]
        for row in par_reader:
            par_rows.append(' '.join(row).split())
    return par_rows

def get_unit(parameter_name:str):
    '''
    #### Input:
        Name of the parameter
    #### Output:
        - One of two
            - Type of unit as an openmm unit
            - Or 1 if the unit is step
    '''
    unit_dict:dict[str] = {
        'FRIC_COEF'         :   unit.picoseconds ** -1,
        'TEMP'              :   unit.kelvin,
        'STEP_SIZE'         :   unit.femtoseconds,
        'PRESSURE'          :   unit.bar,
        'M'                 :   'step',
        'R_MIN'             :   unit.angstrom,
        'R_MAX'             :   unit.angstrom,
        'K_PULL_MIN_EQ'     :   unit.kilocalories_per_mole / unit.angstroms**2,
        'K_PULL_PROD'       :   unit.kilocalories_per_mole / unit.angstroms**2,
        'K_CENTER'          :   unit.kilocalories_per_mole / unit.angstroms**2,
        'DISPLACE_MIN_STEP' :   'step',
        'MIN_STEP'          :   'step',
        'EQ_STEP'           :   'step',
        'SIM_LOOP'          :   'step',
        'REPORT_INTERVAL'   :   'step',
        'SIM_STEP'          :   'step',        
    }
    unit_type = unit_dict[f'{parameter_name.upper()}']
    return unit_type

def get_par(parameter_name:str, main_path:str, get_only_value=False):
    
    paths = get_paths(main_path)
    setup_file = f'{paths["initial_path"]}/setup_file.csv'
    NAME = 0
    VALUE = 1
    UNIT = 2
        
    par_rows = read_parameters(setup_file)
    for row in range(len(par_rows)):
        par = par_rows[row]
        
        if par[NAME] == parameter_name:
            if par[UNIT] == 'step':
                unit_type = 1
                parameter = int(par[VALUE]) * unit_type
            else:
                unit_type = get_unit(parameter_name.upper())
                parameter = float(par[VALUE]) * unit_type
        if get_only_value:
            if par[NAME] == parameter_name:
                if par[UNIT] == 'step':
                    parameter = int(par[VALUE])
                else:
                    parameter = float(par[VALUE])
                
    return parameter    

def get_prompts():
    prompt_dict = {
        'temp'              :   'Please assign the temperature in kelvin:\t',
        'fric_coef'         :   'Please assign the friction coefficient in picoseconds:\t',
        'step_size'         :   'Please assign the step size in femtoseconds:\t',
        'pressure'          :   'Please assing the pressure in bar:\t',
        'M'                 :   'Please assign the number of frames, M, for the umbrella sampling:\t',
        'R_min'             :   'Please assign the lower limit of R0 in angstrom:\t',
        'R_max'             :   'Please assign the upper limit of R0 in angstrom:\t',
        
        'K_pull_min_eq'     :   'Please assign the spring constant, K, for the minimization and equilibration [kcal/mol/angstrom**2]:\t',
        'K_pull_prod'       :   'Please assign the spring constant, K, for the production run [kcal/mol/angstrom**2]:\t',
        'K_center'          :   'Please assign the spring constant, K, for the centering of the host and guest molecules [kcal/mol/angstrom**2]:\t',
        'displace_min_step' :   'Please assign the number of minimization steps after displacement of the host and guest:\t',
        'min_step'          :   'Please assign the number of steps in the mimization process :\t',
        'eq_step'           :   'Please assign the number of steps in the equlibtation process:\t',
        'sim_loop'          :   'Please assign the number of production loops:\t',
        'report_interval'   :   'Please assign how often the distance should be meassured:\t',
        'sim_step'          :   'Please assign how many steps the simulation should run for during each loop:\t',
    }
    return prompt_dict

def get_integer(input_prompt):
    while True:
        num = input(input_prompt)
        try:
            val:int = int(num)
            return val
        except ValueError:
            print('The value must be an integer. Please enter a new number.')

def get_float(input_prompt):
    while True:
        num = input(input_prompt)
        try:
            val:float = float(num)
            return val
        except ValueError:
            print('The value must be an integer. Please enter a new number.')          

def eval_parameters(parameters:list, prompt_dict:dict[str, str]):
    
    parameter_names = [*prompt_dict]
    print_pars = []
    for i in range(len(parameters)):
        unit_type = get_unit(parameter_names[i])
        
        print_pars.append(f'{parameters[i]}\t{unit_type}')
        if parameter_names[i] == 'eq_step':
            eq_steps = parameters[i]
        if parameter_names[i] == 'sim_loop':
            sim_loop = parameters[i]
        if parameter_names[i] == 'sim_step':
            sim_step = parameters[i]
        if parameter_names[i] == 'step_size':
            step_size = parameters[i] * unit_type
    
    eq_run_time = eq_steps * step_size
    prod_run_time = sim_loop * sim_step * step_size
    tot_run_time = eq_run_time + prod_run_time
    
    print(
    f'''
    {"-"*82:82}
    The parameters have been set to the following:
    Integration Parameters:
    Temperature:                                {print_pars[0]}
    Friction coefficient:                       {print_pars[1]}
    Step size:                                  {print_pars[2]}
    Pressure:                                   {print_pars[3]}
    {"-"*82:82}
    
    Umbrella Sampling:
    Number of frames, M:                        {print_pars[4]}
    R_min:                                      {print_pars[5]}
    R_max:                                      {print_pars[6]}
    
    {"-"*82:82}
    
    Spring constants, K:
    K_pull_min_eq:                              {print_pars[7]}
    K_pull_prod:                                {print_pars[8]}
    K_center:                                   {print_pars[9]}
    {"-"*82:82}
    
    Simulation parameters:
    Displacement minimization steps:            {print_pars[10]}
    Minimization steps:                         {print_pars[11]}
    Equlibration steps:                         {print_pars[12]}
    Report Interval:                            {print_pars[13]}
    Production steps:                           {print_pars[14]}
    
    The equlibration will run for:              {eq_run_time}
    The production will run for:                {prod_run_time}
    The total run time of the simulation is:    {tot_run_time}
    {"-"*82:82}
    '''
    )

def write_pars(parameters:list, prompt_dict:dict, setup_path:str):

    parameter_names = [*prompt_dict]
    write_pars = []
    for i in range(len(parameters)):
        unit_type = str(get_unit(parameter_names[i]))
        print(unit_type)
        write_pars.append(f'{parameter_names[i]:20}\t{str(parameters[i]):10}\t{unit_type:32}\n')
        
    with open(setup_path, 'w') as s_file:
        s_file.writelines(write_pars)

def set_pars(main_path:str):
    paths = get_paths(main_path)
    setup_path = f'{paths["initial_path"]}/setup_file.csv'
           
    print(f'Please be ready to give the parameters for the simulation running in {main_path}')
    
    if os.path.exists(setup_path):
        print('\nA setupfile with the following parameters already exists:\n')
        with open(setup_path) as s_file:
            lines = s_file.readlines()
        for i in range(len(lines)):
            print(lines[i].strip('\n'))
    prompt_dict = get_prompts()
    while True:
        
        exists_prompt = 'Do you wish to change the existing parameters? (Y/N)\t>'
        exists_answer = input(exists_prompt)
        while exists_answer.lower() != 'n':    # Spørger om svaret er 'n'
            if exists_answer.lower() == 'y':       # Hvis svaret er 'y' (yes)
                print('Please be ready to put in the parameters:')            # Printer restart prompt
                break   # Breaker dette while loop og går tilbage til det primære!
            else:
                print('Answer not understood. Please only submit Y or N!\n') 
                exists_answer = input(exists_prompt)
                continue    # Fortsætter det sekundære while-loop.
        else:
            exists_answer.lower == 'n'
            break

        temp:float = get_float(prompt_dict['temp'])
        fric_coef:float = get_float(prompt_dict['fric_coef'])
        step_size:float = get_float(prompt_dict['step_size'])
        pressure:float = get_float(prompt_dict['pressure'])
        
        M:int = get_integer(prompt_dict['M'])
        R_min:float = get_float(prompt_dict['R_min'])
        R_max:float = get_float(prompt_dict['R_max'])
        

        K_pull_min_eq:float = get_float(prompt_dict['K_pull_min_eq'])     # kcal/mol/angstrom**2
        K_pull_prod:float = get_float(prompt_dict['K_pull_prod'])         # kcal/mol/angstrom**2
        K_center:float = get_float(prompt_dict['K_center'])               # kcal/mol/angstrom**2

        displace_min_step:int = get_integer(prompt_dict['displace_min_step'])
        min_step:int = get_integer(prompt_dict['min_step'])
        eq_step:int = get_integer(prompt_dict['eq_step'])
        sim_loop:int = get_integer(prompt_dict['sim_loop'])
        report_interval:int = get_integer(prompt_dict['report_interval'])
        sim_step:int = get_integer(prompt_dict['sim_step'])
                
        parameters:list = [temp, fric_coef, step_size, pressure,
                    M, R_min, R_max, 
                    K_pull_min_eq, K_pull_prod, K_center,
                    displace_min_step, min_step, eq_step, sim_loop, report_interval, sim_step]
        
        eval_parameters(parameters, prompt_dict)
        answer:str = input('Do you wish to change anything before writing the paramaters to the setup file? (Y/N)\t>')
        while answer.lower() != 'n':    # Spørger om svaret er 'n'
            if answer.lower() == 'y':       # Hvis svaret er 'y' (yes)
                print('Please submit the parameters again.\n Please take care to do it properly this time')            # Printer restart prompt
                break   # Breaker dette while loop og går tilbage til det primære!
            else:
                print('Answer not understood. Please only submit Y or N!\n') 
                answer = input('Do you wish to change anything? (Y/N)\t>')
                continue    # Fortsætter det sekundære while-loop.
        else:
            write_pars(parameters, prompt_dict, setup_path)  # Det primære while loop er accepteret.
            break

def set_initial_molecules(main_path:str):
    # sub_main = f'{main_path}/{main_path}'
    print(main_path)
    paths = get_paths(main_path)
    molecule_setup= f'{paths["initial_path"]}/initial_molecules.txt'
    print(molecule_setup)
    copy_to_initial_path = paths['initial_path']

    file_extentions = ['SDF', 'SDF', 'PDB']
    molecule_type = ['HOST', 'GUEST', 'COMPLEX']
    header = f'{"Molecule":10}\t{"Initial_file"}\n'

    m_lines = list()
    m_lines.append(header)

    print('The files here are the initial files. They MUST exist in the "in" folder. \nLeave out filename extension.')
    for i in range(3):
        prompt = f'Please enter the name of the {file_extentions[i]} file containing the {molecule_type[i]} molecule:\t'
        
        while True:        
            file = input(prompt)
            in_file = f'in/{file}.{file_extentions[i].lower()}'
            
            if os.path.exists(in_file):
                copy_file = copy(in_file, copy_to_initial_path)
                # write_path = 
                m_lines.append(f'{molecule_type[i].lower():10}\t{copy_file}\n')
                break
            else:
                print('Molecule does not exist. Check spelling')
                continue

    with open(molecule_setup, 'w') as m_files:
        m_files.writelines(m_lines)



def get_ini_molecule_path(molecule:str, main_path:str):
    
    paths = get_paths(main_path)
    molecule_setup = f'{paths["initial_path"]}/initial_molecules.txt'

    with open(molecule_setup, 'r') as m_files:
        lines = m_files.readlines()
    molecules:dict[str,] = {}
    for i in range(len(lines)):
        tokens = lines[i].split()
        if i > 0:
            molecules[f'{tokens[0]}'] = f'{tokens[1]}'    
    return molecules[molecule]

def get_R0(main_path:str, for_shell:bool=False):

    M = get_par('M', main_path)
    R_min = float(str(get_par('R_min', main_path)).split()[0])
    R_max = float(str(get_par('R_max', main_path)).split()[0])
    R0 = np.around(np.linspace(R_min, R_max, M, endpoint=False), 3)

    if for_shell:
        R0 = tuple(R0)
        R0 = str(R0).replace(',', '')
    return R0


def get_R0_index(main_path:str, for_sbatch_array:bool=False, for_bash_loop:bool=False):
    M =  get_par('M', main_path)
    
    R0_index = list(range(M))
    if for_sbatch_array:
        start = R0_index[0]
        stop = R0_index[-1]
        R0_index = f'{start}-{stop}'
    if for_bash_loop:
        start = R0_index[0]
        stop = R0_index[-1]
        R0_index = f'{start} {stop}'
    return R0_index

def copy_py_file(py_file, main_path):
    copy_to_main_path = main_path
    python_script = py_file
    while True:
        in_file = f'{python_script}'
        if os.path.exists(in_file):
            copy_file = copy(in_file, copy_to_main_path)
            print(f'I have copied {python_script:35} {"-->":11} {copy_file:<40}')
            break
        else:
            print('Python script does not exist! Check spelling')
            continue

def get_sub2(main_path:str):
    sub2 = 'C_herc.sh'
    python_script = 'C0.py'
    
    print(f'{python_script} is submitted for {sub2} (creating and setting up the initial files):\t')
    
    R0 = get_R0(main_path, for_shell=True)
    R0_index = get_R0_index(main_path, for_bash_loop=True)
    
    sub2_dict = {
        '0'     :   '#!/usr/bin/env bash',
        '1'     :  f'#SBATCH --job-name C_{main_path}',
        '2'     :   '#SBATCH --partition teach',
        '3'     :   '#SBATCH --time 05:00:00',
        '5'     :   '#SBATCH --ntasks-per-node=1',
        '6'     :   '#SBATCH --cpus-per-task=4',
        '8'     :   '',
        '9'     :   'export PATH=/usr/local/cuda-11.1/bin:$PATH',
        '10'    :   'export LD_LIBRARY_PATH=/usr/local/cuda-11.1/lib64:$LD_LIBRARY_PATH',
        '11'    :   'export OPENMM_CPU_THREADS=$SLURM_CPUS_PER_TASK',
        '12'    :   'source activate pyth39',
        '13'    :   '',
        '14'    :   'ID=$SLURM_ARRAY_TASK_ID',
        '15'    :  f'DISTANCES={R0}',
        '16'    :   '',
        '17'    :  f'for ID in `seq {R0_index}`',
        '18'    :   'do',
        '19'    :  f'    python {python_script} {main_path} ${"{DISTANCES[$ID]}"} $ID',
        '20'    :   'done'
    }

    keys = sub2_dict.keys()
    lines = []
    for key in keys:
        print(f'{sub2_dict[key]}')
        lines.append(f'{sub2_dict[key]}\n')
    
    with open(sub2, 'w') as sub2_file:
        sub2_file.writelines(lines)


def get_sub_AI3(main_path:str):
    
    sub3 = 'D_AIcloud.sh'
    python_script = 'D0_4.py'
    
    print(f'{python_script} is submitted for {sub3} (creating and setting up the initial files):\t')
    
    R0 = get_R0(main_path, for_shell=True)
    R0_index = get_R0_index(main_path, for_sbatch_array=True)
    time_hours = f'{1:02}'
    time_min = f'{30:02}'

    sub3_dict = OrderedDict(
        [
            ('0'     ,   '#!/usr/bin/env bash'),
            ('1'     ,  f'#SBATCH --job-name D_{main_path}'),
            ('2'     ,   '#SBATCH --partition batch'),
            ('3'     ,  f'#SBATCH --time {time_hours}:{time_min}:00'),
            ('4a'    ,   '#SBATCH --gres=gpu:1'),
            ('4b'    ,   '#SBATCH --ntasks-per-node=1'),
            ('4c'    ,   '#SBATCH --cpus-per-task=6'),
            ('4d'    ,   '#SBATCH --qos=allgpus'),
            ('7'     ,  f'#SBATCH --array={R0_index}'),
            ('8'     ,   ''),
            ('9'     ,   'RUN_SH_FILE=run_$SLURM_JOB_NAME.sh'), 
            ('10'    ,   'python3 make_run_sh.py yes $RUN_SH_FILE'),
            ('11'    ,   'chmod +x $RUN_SH_FILE'),
            ('12'    ,   ''),
            ('13'    ,   ''),
            ('14'    ,   'ID=$SLURM_ARRAY_TASK_ID'),
            ('15'    ,  f'DISTANCES={R0}'),
            ('16'    ,   ''),
            ('17'    ,   ''),
            ('18'    ,   ''),
            ('19'    ,  f'MAIN_PATH={main_path}'),
            ('20'    ,   'echo "The main path is $MAIN_PATH. All further simulation will take place here"'),
            ('21'    ,  f'singularity exec ~/py_open_v1.5.img ./$RUN_SH_FILE {python_script} $MAIN_PATH ${"{DISTANCES[$ID]}"} $ID'),
        ]
    )
    
    keys = sub3_dict.keys()
    lines = []
    print('_'*100)
    for key in keys:
        print(f'{sub3_dict[key]}')
        lines.append(f'{sub3_dict[key]}\n')
    print('_'*100, '\n\n'*2)
    
    with open(sub3, 'w') as sub3_file:
        sub3_file.writelines(lines)

def get_sub_herc3(main_path:str):
    
    sub3 = 'D_hercules.sh'
    python_script = 'D0_4.py'
    
    print(f'{python_script} is submitted for {sub3} (creating and setting up the initial files):\t')
    
    R0 = get_R0(main_path, for_shell=True)
    R0_index = get_R0_index(main_path, for_sbatch_array=True)
    time_hours = f'{1:02}'
    time_min = f'{30:02}'

    sub3_dict = OrderedDict(
        [
            ('0'     ,   '#!/usr/bin/env bash'),
            ('1'     ,  f'#SBATCH --job-name D_{main_path}'),
            ('2'     ,   '#SBATCH --partition batch'),
            ('3'     ,  f'#SBATCH --time {time_hours}:{time_min}:00'),
            ('4a'    ,   '#SBATCH --gres=gpu:1'),
            ('4b'    ,   '#SBATCH --ntasks-per-node=1'),
            ('4c'    ,   '#SBATCH --cpus-per-task=6'),
            ('4d'    ,   '#SBATCH --qos=allgpus'),
            ('7'     ,  f'#SBATCH --array={R0_index}'),
            ('8'     ,   ''),
            ('9'     ,   'export PATH=/usr/local/cuda-11.1/bin:$PATH'), 
            ('10'    ,   'export LD_LIBRARY_PATH=/usr/local/cuda-11.1/lib64:$LD_LIBRARY_PATH'),
            ('11'    ,   'export OPENMM_CPU_THREADS=$SLURM_CPUS_PER_TASK'),
            ('12'    ,   ''),
            ('13'    ,   ''),
            ('14'    ,   'SLURM_ID=$SLURM_ARRAY_TASK_ID'),
            ('15'    ,  f'DISTANCES={R0}'),
            ('16'    ,   ''),
            ('17'    ,   ''),
            ('18'    ,   ''),
            ('19'    ,  f'MAIN_PATH={main_path}'),
            ('20'    ,   'echo "The main path is $MAIN_PATH. All further simulation will take place here"'),
            ('21'    ,  f'python {python_script} $MAIN_PATH ${"{DISTANCES[$SLURM_ID]}"} $SLURM_ID'),
        ]
    )
    
    keys = sub3_dict.keys()
    lines = []
    print('_'*100)
    for key in keys:
        print(f'{sub3_dict[key]}')
        lines.append(f'{sub3_dict[key]}\n')
    print('_'*100, '\n\n'*2)
    
    with open(sub3, 'w') as sub3_file:
        sub3_file.writelines(lines)
