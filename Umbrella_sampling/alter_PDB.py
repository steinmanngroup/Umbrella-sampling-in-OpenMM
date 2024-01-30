from sys import stdout

import numpy as np
import matplotlib.pyplot as plt

import csv
from tqdm import tqdm
import rmsd

from openmm.openmm import unit
from openmm.app import PDBFile
from openmm.app import Simulation
from openmm.app import PDBReporter
from openmm.app import StateDataReporter

from create_CCBForce_system import add_centering_force_system
from geometry import dist, normalize_3dvector

#%% Get the lines and its index
def get_pdb_lines(pdb_file:str):
    '''
    Input: 
        PDBFile.
        The file can contain a host, guest and water. '
        It does not take any kind of amino acid and the host must have the lowest indexing number if 
        both at host and guest molecule is present.
    Output:
        The lines of the PDBfile
    '''
    in_pdb_file = pdb_file
    with open(in_pdb_file, 'r') as in_file:
        lines:list[str] = in_file.readlines()
    return lines

def get_atom_index(lines:list[str]):
    '''
    Input:
        A set of lines from a PDBfile.
        Should be combined with the function get_pdb_lines()

    Output:
        The indicies of all lines in the PDBfile starting with ATOM or HETATM
    '''
    atom_index:list[int] = []
    for i, line in enumerate(lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_index.append(i)
    return atom_index

def get_atom_lines_from_pdb(pdb_file:str):
    lines = get_pdb_lines(pdb_file)
    atom_index = get_atom_index(lines)
    
    atom_lines:list[str] = []
    for i in atom_index:
        atom_line = lines[i].split()
        atom_lines.append(atom_line)
    return atom_lines

#%% Use the lines to find the 'tokens' and format these for printing
def get_token_index(token_name:str):
    
    '''
    Input:
        Name of token in PDBfile in string format
        
        Possible inputs:
            atom        In general either HETATM or ATOM
            serial:     gives the number of the atom in the pdb in the specific model
            name:       atomname in the pdb. This is unique for the specific residue
            altLoc:     None
            resName:    Name of the residue
            chainID:    Id of the chain. There can be more of the same residue but the chain should be different
            resSeq:     Sequence of the residue chain
            iCode:      None
            x:          x-coordinate of the atom
            y:          y-coordinate of the atom
            z:          z-coordinate of the atom
            occupancy:  Not sure
            tempFactor: Not sure
            element:    Element i.e. H (hydrogen), C (carbon), O (oxygen) etc.
            charge:     Charge of the atom

    Syntax Source:
        https://pyformat.info/

    Return:
        Index of token
    '''
    pdb_token_dict:dict[str, int] = {
        'atom'      :   0,
        'serial'    :   1,
        'name'      :   2,
        'altLoc'    :   None,
        'resName'   :   3,
        'chainID'   :   4,
        'resSeq'    :   5,
        'iCode'     :   None,
        'x'         :   6,
        'y'         :   7,
        'z'         :   8,
        'occupancy' :   9,
        'tempFactor':   10,
        'element'   :   11,
        'charge'    :   12
        }
    token_index:int = pdb_token_dict[token_name]
    return token_index

def format_tokens(tokens:list[str]):
    '''
    Input:
        The tokens found in each atom-line in a PDBfile.
        A PDBfile is spaced based on the colunmns in the range 1 to 80.
        Each token is specified to hold a specific number of coloulmns according to: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM.

        The coloumns are specified by python formatting: https://pyformat.info/

        format_tokens is combined with the function get_new_lines as this function finds the tokens.
    Output:
        The formatted atoms
    
    '''
    atom = f'{tokens[0]:<6}'                # 1 - 6
    serial = f'{tokens[1]:>5}'              # 7 - 11
    name = f'  {tokens[2]:>3}'              # 13 - 16
    altLoc = f' '                           # 17 (empty due to current pdb format)
    resName = f'{tokens[3]:>1}'             # 18 - 20   
    chainID = f' {tokens[4]:>1}'            # 22 (Character)
    resSeq = f'{tokens[5]:>4}'              # 23 - 26
    iCode = f' '                            # 27 (empty due to current pdb format)
    x = f'   {tokens[6]:>8.3f}'             # 31 - 38, Real(8.3)
    y = f'{tokens[7]:>8.3f}'                # 39 - 46, Real(8.3)
    z = f'{tokens[8]:>8.3f}'                # 47 - 54, Real(8.3)
    occupancy = f'{tokens[9]:>6.2f}'        # 55 - 60, Real(6.2)
    tempFactor = f'{tokens[10]:>6.2f}  '    # 61 - 66, Real(6.2)    
    element = f'{tokens[11]:>10}'           # 77 - 78
    charge = f''                            # 79 - 80 (empty due to current pdb format)

    formatted_tokens:str = f'{atom}{serial}{name}{altLoc}{resName}{chainID}{resSeq}{iCode}{x}{y}{z}{occupancy}{tempFactor}{element}{charge}\n'
    
    return formatted_tokens

def get_residue_name_from_pdb(pdb_file:str):
    '''
    Input:
        Takes a pdb file and finds all residues in it.
        
    Output:
        Returns a list of residue names in the PDBFile.
    '''
    atom_lines = get_atom_lines_from_pdb(pdb_file)
    residue_names_in_pdb = []
    for index, atom_line in enumerate(atom_lines):       
        residue_names_in_pdb.append(atom_line[get_token_index('resName')])
        
    residue_names = []
    for residue_name in residue_names_in_pdb:
        if residue_name not in residue_names:
            residue_names.append(residue_name)
    return residue_names

def get_indicies(pdb_file:str, call_indexing:bool=True, see_indices:bool=False, get_key_names:bool=False):
    lines = get_pdb_lines(pdb_file)
    line_dict:dict[str, int] = {}
    keys:list[str] = []
    for i in range(len(lines)):
        tokens = lines[i].split()
        key:str = tokens[0]
        # print(tokens)
        if key not in line_dict:
            line_dict[key] = [None]
            keys.append(key)
            
        if key == tokens[0]:
            line_dict[key].append(i)
        
        for value_none in line_dict.values():
            if None in value_none:
                value_none.remove(None)
    # print(keys)
    
    atom_index:list[int] = []
    remark_index:list[int] = []
    title_index:list[int] = []
    cryst1_index:list[int] = []
    end_index:list[int] = []
    ter_index:list[int] = []
    conect_index:list[int] = []
    model_index:list[int] = []
    endmdl_index:list[int] = []

    for key in keys:
        for i in line_dict[key]:
            if key == 'TITLE':                      title_index.append(i)
            if key == 'REMARK':                     remark_index.append(i)
            if key == 'CRYST1':                     cryst1_index.append(i)
            if key == 'MODEL':                      model_index.append(i)
            if key == 'HETATM' or key == 'ATOM':    atom_index.append(i)
            if key == 'TER':                        ter_index.append(i)
            if key == 'ENDMDL':                     endmdl_index.append(i)
            if key == 'CONECT':                     conect_index.append(i)
            if key == 'END':                        end_index.append(i)

    ter_host_index = None
    ter_guest_index = None
    ter_solvent_index = None
    if len(ter_index) == 0:
        pass
    elif len(ter_index) > 0:
        ter_host_index:list[int] = [ter_index[0]]
        ter_guest_index:list[int] = [ter_index[1]]
        atom_host_index:list[int] = []
        atom_guest_index:list[int] = []
    
        if len(ter_index) > 2:
            ter_solvent_index:list[int] = [ter_index[2]]
            atom_solvent_index:list[int] = []    
        
    
        for i in atom_index:
            if i < ter_host_index[0]:
                atom_host_index.append(i)
            if i > ter_host_index[0] and i < ter_guest_index[0]:
                atom_guest_index.append(i)
                if len(ter_index) > 2:
                    if i > ter_guest_index[0] and i < ter_solvent_index[0]:
                        atom_solvent_index.append(i)
        
        line_dict['ATOM_HOST_INDEX'] = atom_host_index
        line_dict['TER_HOST_INDEX'] = ter_host_index
        line_dict['ATOM_GUEST_INDEX'] = atom_guest_index
        line_dict['TER_GUEST_INDEX'] = ter_guest_index
        if len(ter_index) > 0:
            insert_keys:list[str] = ['ATOM_HOST_INDEX', 'TER_HOST_INDEX', 'ATOM_GUEST_INDEX', 'TER_GUEST_INDEX']
        if len(ter_index) > 2:
            line_dict['ATOM_SOLVENT_INDEX'] = atom_solvent_index
            line_dict['TER_SOLVENT_INDEX'] = ter_solvent_index
            insert_keys:list[str] = ['ATOM_HOST_INDEX', 'TER_HOST_INDEX', 'ATOM_GUEST_INDEX', 'TER_GUEST_INDEX', 'ATOM_SOLVENT_INDEX', 'TER_SOLVENT_INDEX']
        for i in range(len(keys)):
            # print('finding the keys to remove (HETATM, TER)',keys[i], i, len(keys))
            if keys[i] == 'HETATM' or keys[i] == 'ATOM':
                atom_position = i
            if keys[i] == 'TER':
                ter_position = i
        
        keys.pop(ter_position)
        keys.pop(atom_position)
            
        for key in insert_keys:
            keys.insert(atom_position, key)
            atom_position += 1 
        # print(keys)
        index_list:list[list[int]] = []
        # Title
        if 'TITLE' in keys:
            index_list.append(remark_index)
            
        # Remarks
        if 'REMARK' in keys:
            index_list.append(remark_index)

        # Model
        if 'MODEL' in keys:
            index_list.append(model_index)
        
        # Cryst1
        if 'CRYST1' in keys:
            index_list.append(cryst1_index)

        # Host
        if 'ATOM_HOST_INDEX' in keys:
            index_list.append(atom_host_index)
        if 'TER_HOST_INDEX' in keys:
            index_list.append(ter_host_index)

        # Guest
        if 'ATOM_GUEST_INDEX' in keys:
            index_list.append(atom_guest_index)
        if 'TER_GUEST_INDEX' in keys:
            index_list.append(ter_guest_index)
        
        # Solvent
        if 'ATOM_SOLVENT_INDEX' in keys:
            index_list.append(atom_solvent_index)
        if 'TER_SOLVENT_INDEX' in keys:
            index_list.append(ter_solvent_index)
        
        # Endmdl
        if 'ENDMDL' in keys:
            index_list.append(endmdl_index)
            
        # Conect
        if 'CONECT' in keys:
            index_list.append(conect_index)
            
        #End
        if 'END' in keys:
            index_list.append(end_index)
    # print('length of keys: ', len(keys))
    if call_indexing:
        for i in range(len(keys)):
            
            if see_indices == False:
                if i == 0:
                    print('\nThe PDBFile has the following location and keys')
                    print(f'PDB-location:\t{pdb_file}')
                    print(f'{"_":_>48}')
                    print(f'{"Key index":<9}{"":3}{"Key":20}')
                    print(f'{"_":_>48}')
                print(f'{i:9}{"":3}{keys[i]:20}')
                
                if i == len(keys)-1:
                    print(f'{"_":_>48}\n')
                    print('* Set "see_indices" to True to include indicies in the output.\n')
               
            if see_indices:
                if i == 0:
                    print('\nThe PDBFile has the following location, keys and indices')
                    print(f'PDB-location:\t{pdb_file}')
                    print(f'{"_":_>114}')
                    print(f'{"Key index":<9}{"":3}{"Key":20}{"":3}{"Indices"}')
                    print(f'{"_":_>114}')

                k = 0
                k_step = 15
                while k < len(index_list[i]):
                    indices_to_print = index_list[i][k:  k+k_step]
                    if k < k_step:
                        print(f'{i:9}{"":3}{keys[i]:20}{"|":3}{indices_to_print}')
                    if k >= k_step:
                        print(f'{"":32}{"|":3}{indices_to_print}')
                    k += k_step
                if i == len(keys)-1:
                    print(f'{"_":_>114}\n')
    if get_key_names:
        return index_list, keys
    return index_list

def get_keys_index(keys:list[str]):
    '''
    Takes a list of keys and produces an index of these.
    Specifically it is setup to find the ATOM_HOST_INDEX and ATOM_GUEST_INDEX
    '''
    key_index = list(range(2))
    for i in range(len(keys)):
        if keys[i] == 'ATOM_HOST_INDEX':
            key_index[0] = i
        if keys[i] == 'ATOM_GUEST_INDEX':
            key_index[1] = i
    return key_index

def get_lines_from_index(index:list[int], lines:list[str]):
    '''
    Input:
        Index and the lines to be written
    
    Example:
        index_list = get_indecies(pdb_file, see_indices=False)\n
        lines = get_pdb_lines(pdb_file)
        
        Write the first index of index_list:
        with open('file', 'w') as file:
            file.writelines(write_lines_to_file(index_list[0], lines))
    '''
    lines_from_index:list[str] = []
    for i in index:
        lines_from_index.append(lines[i])
    return lines_from_index

def get_host_info(pdb_file:str, token_name:str):
    token_name_list_str:list[str] = ['atom', 'name', 'element', 'chainID']
    token_name_list_float:list[str] = ['x', 'y', 'z', 'charge', 'tempFactor', 'occupancy']
    
    if token_name not in token_name_list_float + token_name_list_str + ['serial' ,'resSeq']:
        raise Exception(f'{token_name} is not a viable token name!')
    else:
        host_info:list = []
        lines:list[str] = get_pdb_lines(pdb_file)
        index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
        key_index = get_keys_index(keys)[0]
        host_index = index_list[key_index]
        host_lines = get_lines_from_index(host_index, lines)
        

        if token_name in token_name_list_str:    
            for i in range(len(host_lines)):
                host_info.append(host_lines[i].split()[get_token_index(token_name)])
        
        if token_name in token_name_list_float:
            for i in range(len(host_lines)):
                host_info.append(float(host_lines[i].split()[get_token_index(token_name)]))
            
        if token_name == 'serial':
            for i in range(len(host_lines)):
                host_info.append(int(host_lines[i].split()[get_token_index(token_name)]) - host_index[0] + 1)
                
        if token_name == 'resSeq':
            for i in range(len(host_lines)):
                host_info.append(int(host_lines[i].split()[get_token_index(token_name)]))
        return host_info


def get_guest_info(pdb_file:str, token_name:str):
    token_name_list_str:list[str] = ['atom', 'name', 'element', 'chainID']
    token_name_list_float:list[str] = ['x', 'y', 'z', 'charge', 'tempFactor', 'occupancy']
    
    if token_name not in token_name_list_float + token_name_list_str + ['serial' ,'resSeq']:
        raise Exception(f'{token_name} is not a viable token name!')
    else:
        guest_info:list = []
        lines:list[str] = get_pdb_lines(pdb_file)
        index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
        key_index = get_keys_index(keys)
        host_index = index_list[key_index[0]]
        
        guest_index = index_list[key_index[1]]
        guest_lines = get_lines_from_index(guest_index, lines)
        
        if token_name in token_name_list_str:    
            for i in range(len(guest_lines)):
                guest_info.append(guest_lines[i].split()[get_token_index(token_name)])
        
        if token_name in token_name_list_float:
            for i in range(len(guest_lines)):
                guest_info.append(float(guest_lines[i].split()[get_token_index(token_name)]))
            
        if token_name == 'serial':
            for i in range(len(guest_lines)):
                guest_info.append(int(guest_lines[i].split()[get_token_index(token_name)]) - host_index[0])
                
        if token_name == 'resSeq':
            for i in range(len(guest_lines)):
                guest_info.append(int(guest_lines[i].split()[get_token_index(token_name)]))
        return guest_info

def get_host_serial(pdb_file:str):
    host_serial:list[int] = []
    
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
    key_index = get_keys_index(keys)[0]
    host_index = index_list[key_index]
    host_lines = get_lines_from_index(host_index, lines)
    
    for i in range(len(host_lines)):
        host_serial.append(int(host_lines[i].split()[get_token_index('serial')]) - host_index[0] + 1)
    return host_serial

def get_guest_serial(pdb_file:str):
    guest_serial:list[int] = []
    
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
    key_index = get_keys_index(keys)
    host_index = index_list[key_index[0]]
    
    guest_index = index_list[key_index[1]]
    guest_lines = get_lines_from_index(guest_index, lines)
    for i in range(len(guest_lines)):
        guest_serial.append(int(guest_lines[i].split()[get_token_index('serial')]) - host_index[0])
    return guest_serial

def get_RC(pdb_file):
    '''
    C3 and C5 from the sugar molecules are used to find P_{top} and P_{bottom} and the center of mass of these two.
    RC is the vector describing the direct way through the two center of mass and the line of the Reaction Coordinate (RC)
    
    Input:
        Index list of atoms from PDBFile (get_indicies)
        Lines of PDBFile (get_pdb_lines)

    Output:
        The normalized RC:
        RC_norm = [x, y, z]
    '''
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
       
    atom_host_lines = []
    key_index = get_keys_index(keys)[0]

    for i in index_list[key_index]:
        atom_host_lines.append(lines[i])
        
    
    atoms = ['C3', 'C5']
    P = list(range(3))
    RC_diff = list(range(len(atoms)))
    # print('RC_diff: ', RC_diff)
    # print('atoms: ', atoms)
    for j, atom in enumerate(atoms):
        atom_types = []
        coordinates = []
        # print('outside',j, atom)
        for k in range(len(atom_host_lines)):
            tokens = atom_host_lines[k].split()
            
            if tokens[2] == atom:
                # print(f'index {k:03d} and accepted token:{"":3} {str(tokens)}')
                # print('inside', j, atom)
                x = float(tokens[get_token_index('x')])
                y = float(tokens[get_token_index('y')])
                z = float(tokens[get_token_index('z')])
                coordinates.append([x, y, z])
                atom_types.append(tokens[get_token_index('element')])
       
        P = rmsd.get_cm(atom_types, coordinates)
        RC_diff[j] = P
        
    RC:list = RC_diff[0] - RC_diff[1]
    RC_norm = np.array(normalize_3dvector(RC[0], RC[1], RC[2]))
    
    print(f'\nRC:{"":11} [{RC[0]:>03.5f}, {RC[1]:>02.5f}, {RC[2]:>02.5f}]')
    print(f'Normalized RC: [{RC_norm[0]:>03.5f}, {RC_norm[1]:>02.5f},  {RC_norm[2]:>03.5f}]\n')
    
    return RC_norm

def plot_in_3d(points:np.matrix, title:str):
    x = []
    y = []
    z = []
    for i in range(len(points)):
        x.append(points[i][0])
        y.append(points[i][1])
        z.append(points[i][2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(x, y, z, cmap='Reds')
    ax.view_init(0, 0)
    ax.text2D(0.05, 0.95, f'{title}', transform=ax.transAxes)
    fig.savefig(title + '.png')


def move_centered_complex_to_orego(lines, index_list, key_index):
    '''
    Moving host and guest molecule to orego by substracting the center of mass of the host molecule.
    
    Input:
        Lines and index of lines
        
    Output:
        Lines with the new coordinates.
    
    '''

    for j in key_index:
        atom_lines:list = []
        atom_types = []
        P = []

        for i in index_list[j]: 
            atom_lines.append(lines[i])
        
        for i in range(len(atom_lines)):
            tokens = atom_lines[i].split()
            
            x = float(tokens[get_token_index('x')])
            y = float(tokens[get_token_index('y')])
            z = float(tokens[get_token_index('z')])
            P.append([x, y, z])
            atom_types.append(tokens[get_token_index('element')])
        
        if j == key_index[0]:
            center_mass = rmsd.get_cm(atom_types, P)
            new_P_host = P - center_mass
            
        if j == key_index[1]:
            new_P_guest = P - center_mass
            
    return new_P_host, new_P_guest

def get_centered_complex_lines_at_orego(pdb_file):
    
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)
    key_index = get_keys_index(keys)
    new_coordinates = move_centered_complex_to_orego(lines, index_list, key_index)
    
    k = 0
    for j in key_index:
        atom_lines = []
        new_lines = []
        for i in index_list[j]: 
            atom_lines.append(lines[i])
        for i in range(len(atom_lines)):
            tokens = atom_lines[i].split()
        
            tokens[get_token_index('occupancy')] = float(tokens[get_token_index('occupancy')])
            tokens[get_token_index('tempFactor')] = float(tokens[get_token_index('tempFactor')])
            
            tokens[get_token_index('x')] = new_coordinates[k][i][0]
            tokens[get_token_index('y')] = new_coordinates[k][i][1]
            tokens[get_token_index('z')] = new_coordinates[k][i][2]
            
            formatted_tokens = format_tokens(tokens)
            new_lines.append(formatted_tokens)
        k += 1
        if j == key_index[0]: host_new_lines = new_lines
        if j == key_index[1]: guest_new_lines = new_lines
    return host_new_lines, guest_new_lines

def write_centered_complex_to_orego(pdb_file:str, out_path:str):
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=True, see_indices=True, get_key_names=True)

    pdb_at_orego = f'{out_path}/orego.pdb'
    with open(pdb_at_orego, 'w') as file:
        for i in range(len(keys)):        
            if keys[i] == 'ATOM_HOST_INDEX':
                file.writelines(get_centered_complex_lines_at_orego(pdb_file)[0])
            if keys[i] == 'ATOM_GUEST_INDEX':
                file.writelines(get_centered_complex_lines_at_orego(pdb_file)[1])
            if keys[i] not in ['ATOM_HOST_INDEX', 'ATOM_GUEST_INDEX']:
                file.writelines(get_lines_from_index(index_list[i], lines))
    print(f'\nPDBFile with centered host and guest molecules are made: \n')
    print(f'\tLocation of the initial PDBFile: {pdb_file}')
    print(f'\tLocation of the new PDBFile:     {pdb_at_orego}\n\n')
    return pdb_at_orego

def displace_guest(pdb_file:str, R0:list[float], R0_index, displace_path:str):
    '''
    The guest molecule is displaced by a distance defined by R0:
        R0 :list[float] = np.linspace(start=R_min, stop=R_max, num=M, endpoint=False)
    Each step displaces the guest relative to the starting point, i.e. the Reaction Coordinate (RC).
    For each step a new PDBFile is written, representing the window in an Umbrella Sampling (US).
    
    Input:
        PDBFile containing host and guest molecule.
        R0 is a list of floats increasing step by step.
        Outpath should point to the mainfolder
        
    Return:
        Path to the PDBFile with the displaced guest molecule
    '''
    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, get_key_names=True)

    RC_norm = get_RC(pdb_file)
    
    guest_key_index = get_keys_index(keys)[1]
    atom_lines:list = []
    for i in index_list[guest_key_index]: 
        atom_lines.append(lines[i])

    window_path = f'{displace_path}'
       
    step_length = RC_norm * R0
    new_lines = []
    for i in range(len(atom_lines)):
        tokens = atom_lines[i].split()
        x = float(tokens[get_token_index('x')])
        y = float(tokens[get_token_index('y')])
        z = float(tokens[get_token_index('z')])
        P = step_length + [x, y, z]
        
        tokens[get_token_index('x')] = P[0]
        tokens[get_token_index('y')] = P[1]
        tokens[get_token_index('z')] = P[2]
        tokens[get_token_index('occupancy')] = float(tokens[get_token_index('occupancy')])
        tokens[get_token_index('tempFactor')] = float(tokens[get_token_index('tempFactor')])
    
        formatted_tokens = format_tokens(tokens)
        new_lines.append(formatted_tokens)
    
    pdb_window = f'{R0_index:03d}.pdb'
    pdb_window_path = f'{window_path}/{pdb_window}'
    with open(pdb_window_path, 'w') as file:
        for i in range(len(keys)):        
            if keys[i] == 'ATOM_GUEST_INDEX':
                file.writelines(new_lines)
            if keys[i] not in 'ATOM_GUEST_INDEX':
                file.writelines(get_lines_from_index(index_list[i], lines))        
    print(f'Window {R0_index:03d} is located here:{"":3}{pdb_window_path}{"":3} Displacement factor: R0 = {R0:02.3f}')
    return pdb_window_path
     
def get_distance(pdb_file):
    '''
    Returns the distance between a host and a guest molecule in a pdb.
    Any water in the pdb will be ignored.
    
    Variable:
        pdb_file: 
    '''

    lines = get_pdb_lines(pdb_file)
    index_list, keys = get_indicies(pdb_file, call_indexing=False, see_indices=False, get_key_names=(True))
    keys_index = get_keys_index(keys)


    cm_m:list[float] = list(range(2))
    k = 0
    
    for j in keys_index:
        atom_lines:list[str] = []
        atom_types:list[str] = []
        P = []
        
        for i in index_list[j]: 
            atom_lines.append(lines[i])
        
        for i in range(len(atom_lines)):
            tokens = atom_lines[i].split()
            
            x = float(tokens[get_token_index('x')])
            y = float(tokens[get_token_index('y')])
            z = float(tokens[get_token_index('z')])
            P.append([x, y, z])
            atom_types.append(tokens[get_token_index('element')])
        
        # print(rmsd.get_cm(atom_types, P), k)
        cm_m[k] = (rmsd.get_cm(atom_types, P))
        # print('after cm: ', cm_m)
        k =+ 1
        
    v1, v2 = cm_m
    host_guest_distance:float = dist(v1, v2)
    return host_guest_distance

#%% Running the simulations outside the main file
def get_TER_pdb(initial_pdb, out_path, simulation_parameters):
    '''
    Takes a 'fresh' PDBFile and adds TER lines to it in a format openmm prefers.

    The simulation parameters are a list containing the modeller, system and integrator and must be set up as follows
        simulation_parameters = [modeller, system, integrator]
        
    Returns:
        PDBFile containing the TER lines, which openmm prefers
    '''
    modeller, system, integrator = simulation_parameters
    
    print(f'\nThe integrator is setup with: \n\tTemperature: {integrator.getTemperature()} \t Friction: {integrator.getFriction()} \tStepsize: {integrator.getStepSize()}')

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    ter_pdb = f'{out_path}/ter.pdb'
    simulation.reporters.append(PDBReporter(ter_pdb, reportInterval=1))
    simulation.step(1)
    
    print(f'\nPDBFile with TER (terminating lines) are made: \n')
    print(f'\tLocation of the initial PDBFile: {initial_pdb}')
    print(f'\tLocation of the new PDBFile:     {ter_pdb}\n\n')
    
    return ter_pdb

def get_centered_pdb(pdb_file_ter:str, out_path:str, simulation_parameters:list):
    '''
    simulation parameters are a list containing the modeller, system and integrator and must be set up as follows
        simulation_parameters = [modeller, system, integrator]
        
    Returns:
        PDBFile containing the guest centered inside the host
    '''
    modeller, system, integrator = simulation_parameters
    
    host_serial = get_host_serial(pdb_file_ter)
    guest_serial = get_guest_serial(pdb_file_ter)

    center_force_system = add_centering_force_system(host_serial, guest_serial, system)
    simulation = Simulation(modeller.topology, center_force_system, integrator)
    simulation.context.setPositions(modeller.positions)

    for i in tqdm(range(200)):
        simulation.minimizeEnergy(maxIterations=20)
        
    pdb_file_center = f'{out_path}/center.pdb'
    
    report_steps = 100000
    simulation.reporters.append(PDBReporter(pdb_file_center, reportInterval=report_steps))
    simulation.reporters.append(StateDataReporter(stdout, reportInterval=report_steps/10, step=True,totalSteps=report_steps, progress=True, separator='\t', elapsedTime=True))
    
    K_center = 10000.0 * unit.kilocalories_per_mole / 1.0  * unit.angstroms ** 2
    simulation.context.setParameter('K_center', K_center)
    simulation.step(report_steps)
    
    return pdb_file_center
    
def get_frame0(pdb, production_pdbs_path):
    '''
    Writes the first frame of a pdb file as a new pdb with a single frame in.
    
    Variable:
        Takes <openmm.app.pdbfile.PDBFile> object
        
    Return:
        pdb file with the name 'frame0.pdb'
    '''
    frame0_handle = open(f'{production_pdbs_path}/frame0.pdb', 'w')
    PDBFile.writeFile(topology=pdb.getTopology(), 
                      positions=pdb.getPositions(frame=0), 
                      file=frame0_handle)
    frame0_handle.close()
    return f'{production_pdbs_path}/frame0.pdb'

def get_frames_distances(pdb_file:str):
    pdb = PDBFile(pdb_file)
    frames = pdb.getNumFrames()
    
    frame0 = get_frame0(pdb)
    host_serial = get_host_info(frame0, 'serial')
    guest_serial = get_guest_info(frame0, 'serial')
    host_elements = get_host_info(frame0, 'element')
    guest_elements = get_guest_info(frame0, 'element')

    host_guest_distance:list[float] = list(range(frames))
    
    if frames >= 20:        # If more than 20 frames prints every 5% of the progress
        progres = range(0, frames, int(frames*0.05))
    if frames < 20:        # If less than 20 frames prints all progression
        progres = range(0, frames)
    for i in range(frames):
        if i in progres:
            print(f'Progresion: {i:05d} distances of {frames} ({i/frames*100:05.2f}%)')
            
        pos = pdb.getPositions(frame=i, asNumpy=True)
        host_pos, guest_pos = [], []
        for j in host_serial + guest_serial:
            pos_line = pos[j].in_units_of(unit.angstroms)
            if j in host_serial: 
                host_pos.append(pos_line)
                if j == host_serial[-1]:
                    v1 = rmsd.get_cm(host_elements, host_pos)
            if j in guest_serial:
                guest_pos.append(pos_line)
                if j == guest_serial[-1]:
                    v2 = rmsd.get_cm(guest_elements, guest_pos)
                    host_guest_distance[i] = dist(v1, v2)
    return host_guest_distance

def get_frames_distances2(pdb, serials:list[int], elements:list[int], frames:int):
    print(type(serials[0]), type(serials[0][0]))
    frames = int(frames)
    host_serial, guest_serial = serials
    host_elements, guest_elements = elements

    host_guest_distance:list[float] = list(range(frames))
    
    # if frames >= 20:        # If more than 20 frames prints every 5% of the progress
    #     progres = range(0, frames, int(frames*0.05))
    # if frames < 20:        # If less than 20 frames prints all progression
    #     progres = range(0, frames)
    for i in range(frames):
        # if i in progres:
        #     print(f'Progresion: {i:05d} distances of {frames} ({i/frames*100:05.2f}%)')
            
        pos = pdb.getPositions(frame=i, asNumpy=True)
        host_pos, guest_pos = [], []
        for j in host_serial + guest_serial:
            pos_line = pos[j].in_units_of(unit.angstroms)
            if j in host_serial: 
                host_pos.append(pos_line)
                if j == host_serial[-1]:
                    v1 = rmsd.get_cm(host_elements, host_pos)
            if j in guest_serial:
                guest_pos.append(pos_line)
                if j == guest_serial[-1]:
                    v2 = rmsd.get_cm(guest_elements, guest_pos)
                    host_guest_distance[i] = dist(v1, v2)
    return host_guest_distance


def get_and_write_distances(production_pdb, sim_loop, dist_file_path):
    print(f'Reading the distances of {production_pdb} now')
    progress = range(0, sim_loop, int(sim_loop*0.05))
    res_name = []
    with open(production_pdb, 'r') as n_file:
        res_counter = 0
        for line in n_file:
            if line.startswith('TER'):
                tokens = line.split()
                res_name.append(tokens[2])
                res_counter += 1
            if res_counter == 2:
                break

    with open(dist_file_path, 'w') as dist_file:
        write = csv.writer(dist_file, delimiter='\n')
        with open(production_pdb, 'r') as p_file:
            host_coords = []
            guest_coords = []
            host_elements = []
            guest_elements = []
            j = 0
            k = 0
            count_frames = 0
            cm_m:list[float] = list(range(2))
            for line in p_file:
                
                k += 1
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    tokens = line.split()
                    if tokens[get_token_index('resName')] == res_name[0]:
                        x = float(tokens[get_token_index('x')])
                        y = float(tokens[get_token_index('y')])
                        z = float(tokens[get_token_index('z')])
                        host_elements.append(tokens[get_token_index('element')])
                        host_coords.append([x, y, z])
                        cm_m[0] = (rmsd.get_cm(host_elements, host_coords))
                        # print(host_coords)
                    if tokens[get_token_index('resName')] == res_name[1]:
                        x = float(tokens[get_token_index('x')])
                        y = float(tokens[get_token_index('y')])
                        z = float(tokens[get_token_index('z')])
                        guest_coords.append([x, y, z])
                        guest_elements.append(tokens[get_token_index('element')])
                        cm_m[1] = (rmsd.get_cm(guest_elements, guest_coords))
                        # print(guest_coords)
                if line.startswith('TER'):
                    # print(f'k: {k}')
                    j += 1
                    if j == 2:
                        host_coords.clear()
                        host_elements.clear()
                        guest_coords.clear()
                        guest_elements.clear()
                        v1, v2 = cm_m
                        host_guest_distance = [(dist(v1, v2))]
                        
                        write.writerow(host_guest_distance)
                        
                    if j == 3:
                        count_frames += 1
                        j = 0
                        if count_frames in progress:
                            print(f'Progresion: {count_frames:05d} distances of {sim_loop} ({count_frames/sim_loop*100:03.0f}%)')
