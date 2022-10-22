from sys import argv


from setup import set_initial_molecules, set_pars, set_paths
# from setup import get_sub2, get_sub_AI3

script, main_path = argv

set_paths(main_path=main_path)
set_pars(main_path=main_path)
set_initial_molecules(main_path)

# Setting up sh files
# get_sub2(main_path=main_path)
# get_sub_AI3(main_path=main_path)
