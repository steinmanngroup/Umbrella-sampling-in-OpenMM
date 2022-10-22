"""
Module to create, write and read CustomCentroidBondForce systemfiles.
"""

from pathlib import Path

from openmm.openmm import XmlSerializer
from openmm.openmm import unit
from openmm.openmm import VerletIntegrator
from openmm.openmm import CustomCentroidBondForce

from openmm.app import Simulation
from openmm.app import ForceField
from openmm.app import Modeller
from openmm.app import PDBFile
from openmm.app import NoCutoff

from openmm.vec3 import Vec3

# OpenFF
from openff.toolkit.topology import Molecule

# From SMIRNOFF
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

def get_SMIRNOFF_molecule_forcefield(in_sdf_file):
    '''
    Generates a force field with 
    '''
    smirnoff = SMIRNOFFTemplateGenerator(forcefield='openff-2.0.0.offxml')
    
    # Adding monovalent ions to neutralize any charge
    ion_smiles = ['[Na+]', '[Cl-]']
    for ion_smile in ion_smiles:
        ion = Molecule.from_smiles(ion_smile)
        smirnoff.add_molecules(ion)
    
    offmol_sdf = in_sdf_file
    for i in range(2):
        offmol = Molecule.from_file(offmol_sdf[i])
        smirnoff.add_molecules(offmol)
        
    forcefield = ForceField('amber14-all.xml', 'tip3p.xml', 'amber/GLYCAM_06j-1.xml')
    forcefield.registerTemplateGenerator(smirnoff.generator)
    return forcefield

def smiles_to_sdf(smiles_string:str, out_sdf:str, molecule_name:str):
    '''
    Make a sdf file from a smiles string. The hydrogen cannot be explicit!
    
        smiles_string:  A string representation of a molecule. https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
        
        out_sdf:        Path to the sdf file. Must be a string
        molecule:       Name of the molecule. Must be a string
    '''
    smi_str = smiles_string
    mol_smi = Molecule.from_smiles(smi_str, hydrogens_are_explicit=False, allow_undefined_stereo=True)
    mol_smi.name = f'{molecule_name}'
    mol_smi.visualize()
    
    if '*.sdf' not in out_sdf:
        file_path = out_sdf + '.sdf'
        print(file_path)
    
    Molecule.to_file(mol_smi, file_path=file_path, file_format='SDF')

def add_pull_force_system(host_serial, guest_serial, system):
    pull_force = CustomCentroidBondForce(2, '0.5*K_pull*dR^2; dR=(R-R0); R=distance(g1, g2)')
    pull_force.addGlobalParameter('K_pull', 0.0 * unit.kilocalories_per_mole / 1.0  * unit.angstroms ** 2)
    pull_force.addGlobalParameter('R0', 0.0)
    host = pull_force.addGroup(host_serial)
    guest = pull_force.addGroup(guest_serial)
    pull_force.addBond([host, guest])
    # print(pull_force.getEnergyFunction())
    system.addForce(pull_force)
    return system



def add_centering_force_system(host_serial, guest_serial, system):
    set_center_force = CustomCentroidBondForce(2, 'K_center * distance(g1, g2)')
    set_center_force.addGlobalParameter('K_center', 0.0 * unit.kilocalories_per_mole / 1.0  * unit.angstroms ** 2)
    host = set_center_force.addGroup(host_serial)
    guest = set_center_force.addGroup(guest_serial)
    set_center_force.addBond([host, guest])
    # print(f'{"Energy function of centering force:":}', set_center_force.getEnergyFunction())
    system.addForce(set_center_force)
    return system

def create_modeller(in_pdb_file:str, forcefield, water:bool=False, padding:unit=0.5*unit.nanometers):
    pdb = PDBFile(in_pdb_file)

    modeller = Modeller(pdb.topology, pdb.positions)
    if water:
        modeller.addSolvent(forcefield, padding=padding, model='tip3p')
    return modeller

def create_modeller_box(in_pdb_file, forcefield, water=False, box_side=5.0):
    '''
    Create a modeller with a pdb. 
    
    Optional:
        Adds water a and a boxsize of 5 nm^3 as a standard. Change the value 'box_side' to get a different dimesion.
    '''
    pdb = PDBFile(in_pdb_file)

    modeller = Modeller(pdb.topology, pdb.positions)
    a_vec = Vec3(box_side, 0, 0)
    b_vec = Vec3(0, box_side, 0)
    c_vec = Vec3(0, 0, box_side)
    if water:
        modeller.addSolvent(forcefield, model='tip3p', boxVectors=(a_vec, b_vec, c_vec))
    return modeller


def create_centroid_system(in_pdb_file:str, in_sdf_file:list[str], water=True, padding=0.50*unit.nanometers):
    """
    Creates a system with CustomCentroidBondForce (CCBForce) between two groups.

    Input:
        in_pdb_file: A single PDB-file containing all residues needed for the system.
        in_sdf_file: A list of SDF-files containing two the two structures which needs to be pulled apart.
        water: default is False. Adds tip3p water to the system
    
    Returns:
        The system.
    Verification:
        A simulation-object will be made, before writing the system.xml to verify the syntax energy string.
        If the syntax is improper, the error will be returned.
    """

    forcefield = get_SMIRNOFF_molecule_forcefield(in_sdf_file)
    modeller = create_modeller(in_pdb_file, forcefield, water=water, padding=padding)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)
    
    print('Number of particles: ', system.getNumParticles())
    
    # Should be set up to depend on the residues and possible the chain 
    # instead of finding the interval manually.
    index_host = [ index for index in range(147) ]
    index_guest = [ index for index in range(147, 159)]
    
    pull_force = CustomCentroidBondForce(2, '0.5*K_pull*dR^2; dR=(R-R0); R=distance(g1, g2)')
    pull_force.addGlobalParameter('K_pull', 0.0 * unit.kilocalories_per_mole / 1.0  * unit.angstroms ** 2)
    pull_force.addGlobalParameter('R0', 0.0 * unit.angstroms)
    host = pull_force.addGroup(index_host)
    guest = pull_force.addGroup(index_guest)
    pull_force.addBond([host, guest])
    print(pull_force.getEnergyFunction())
    system.addForce(pull_force)


    set_center_force = CustomCentroidBondForce(2, 'K_center * distance(g1, g2)')
    set_center_force.addGlobalParameter('K_center', 0.0 * unit.kilocalories_per_mole / 1.0  * unit.angstroms ** 2)
    host = set_center_force.addGroup(index_host)
    guest = set_center_force.addGroup(index_guest)
    set_center_force.addBond([host, guest])
    print(set_center_force.getEnergyFunction())
    system.addForce(set_center_force)
    try: # Jeg vil gerne have en kontrol af energifunktionerne, som måske ikke kræver en simulering
        Simulation(modeller.topology, system, integrator=VerletIntegrator(1*unit.femtoseconds))
                
    except ValueError:
        print('An error occured') # Kan det ikke bare være standard error meldingen i stedet?
        raise
    except TypeError:
        print('noget')
        raise
    else:
        print('The system was created succesfully')
        return system


def write_system_to_file(system_variable, write_to_folder:str):
    """
    Writes a system to an xml-file and prints the path
    """
    Path(f'{write_to_folder}').mkdir(parents=True, exist_ok=True)
    system_file_path = f'{write_to_folder}system.xml'
    print(f'Writes the system to an xml-file with the name: {system_file_path}')
    with open(system_file_path, 'w') as file_handle:
        file_handle.write(XmlSerializer.serialize(system_variable))
        print(system_file_path)
    return system_file_path

    

def read_system_from_file(system_file_xml:str):
    """
    Reads an openmm-system from an xml-file and returns it.
    """
    
    system_file_path= system_file_xml
    print('Reads a system to an xml-file with the name: ', {system_file_path})
    with open(system_file_path, 'r') as file_handle:
        xml = file_handle.read()
    system = XmlSerializer.deserialize(xml)
    
    return system


def get_custom_system_forces(system_file_xml:str):
    search_force_energy = '<Force energy='
    search_parameter = '<Parameter default='

    file_name = system_file_xml
    with open(file_name) as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines, start=0):
              
            if search_force_energy in line:
                force_energy_line = i
                print(lines[force_energy_line])
            if search_parameter in line:
                parameter_line = i
                print(lines[parameter_line])

    