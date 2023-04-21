# Umbrella sampling in OpenMM

The purpose of this project was to simulate the guest-host interactions between beta-cyclodextrins (bCD) and a ligand.

## Data
The data sampled in this project is available for download at https://zenodo.org/record/7852102#.ZEJeIHZBwuU.


## Requirements

**Molecular dynamic simulation:**  
OpenMM: `conda install -c conda-forge openmm` [Source](http://docs.openmm.org/latest/userguide/)

Openforcefield `conda install -c conda-forge openmmforcefields` [Source](https://github.com/openmm/openmmforcefields)

**Data treatment:**  
FastMBAR: `pip install -U FastMBAR` [Source](https://fastmbar.readthedocs.io/en/latest/index.html)  
It is recommended to follow their guide to ensure a proper installation. 


## Molecular dynamic simulation
### How to
The directory "in" contains an assorment of sdf and pdb files which should be paired up by name.

`A0.py` creates the initial files and directories and takes `main_path` as its single argument which must be constant for all the following python files (`B0, C0, D0_4`)
`A0.py` copies the the nessary files from the ´in´-directory and sets up the paths at the same time.
When promted either set the parameters or use the ones supplied in the benzene directory (`setup_file.csv`). 


`B0.py` setup the initial pdb-files. It takes the single input `main_path `. The program ensures that the host and guest are centered on top of each other and the starting point for all future simulations is orego of the system. It should be suficcient to run `B0.py` only once per `main_path`

`C0.py` takes three arguments: `main_path R0 R0_index`. R0 is the distance between the host and ligand, while R0_index is the index of the distance: R0 = [0, 0.25, 0.50, 0.75, 1.00, n]; R0_index = [0, 1, 2, 3, 4, n_i]. The Umbrella Sampling is setup here with the distances of `R0` as the molecules are setup with their initial distances.

`D0_4.py` takes the same arguments as `C0.py` `main_path R0 R0_index`. Here the host and guest is solvated in water and the production run is made. At each distance, set by `R0`, a `CustomCentroidForce` is used to keep the distance of `R0`. The distance is controlled by the harmonic potential $E_{pot} = \frac{1}{2} \cdot K \cdot (R-R_0)^2$, where K is the spring constant (usually set to 100 $kJ M^{-1}$ and $R$ is the measured distance.

The measured distances are written to a csv-file in `main_path/production/distance/`.


### SLURM MD-simulation
Run the pythonfile `A0.py`. This will produce two additional files: `C_herc.sh` and `D_hercules.sh`.
Both of these are setup for a specific slurm system but in essence, they automate much of the procedures.
`B0.py` still has to be run seperatly.




## Data treatment 
This specific approach is made for data sampled through umbrella-sampling and is heavily inspired by the FastMBAR [documentation](https://fastmbar.readthedocs.io/en/latest/examples.html). 

### How to
`F0_3.py` must be placed in the same directory as the other scripts. `F0_3.py` takes four arguments: `main_path run_on_gpu high L_times`. `main_path` must be the same as for `D0, C0`, and `D0_4`. `run_on_gpy` can be either 'yes' or 'no', `high` ask at what side of the bCD the ligand is placed, and `L_times` asks how many windows the umbrella sampling should be split into.

With the dataset suggested here, the argument `high` should be set to 'yes', as it is positive.

### SLURM data treatment
In SLURM, run the file `F0_3` and follow the same procedure as for `F0_3.py`.
