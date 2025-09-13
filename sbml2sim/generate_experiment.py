import sys

import s2s
from pathlib import Path
import json
from shutil import copy2

def main():
    sbml_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    sbml = s2s.SBMLDoc(sbml_path)
    
    constants = sbml.get_kinetic_constants()
    
    config_path = Path(output_dir) / "config.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with config_path.open("w") as f:
        f.write("files:\n")
        f.write(f"    model: ./{Path(sbml_path).name}\n")
        f.write(f"    parameters: ./parameters.tsv\n")
        f.write(f"    parameters_costrains: ./constrains.tsv\n")
        f.write("""horizon: 196
absolute_tolerance: 1.0e-12
relative_tolerance: 1.0e-06
time_step: 1

loss_functions:
    stability:
        scale: 10
    transitorial:
        scale: 1

simulator: RoadRunner
solver: nevergrad
random_seed: 104
solver_parameters:
    optimizer: DE
    budget: 100 """)
    
    parameters_file = Path(output_dir) / "parameters.tsv"
    with parameters_file.open("w") as f:
        f.write("parameterId	parameterName	parameterScale	lowerBound	upperBound	nominalValue	estimate	distribution	mu	sigma	granularity\n")
        for constant in constants:
            lower_bound = -20.0 # TODO paremetrize
            upper_bound = 20.0
            granularity = 0.0
            f.write(f"{constant}	{constant}_name	log	{lower_bound}	{upper_bound}	0.0	1	normal	0.0	0.0	{granularity}\n")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    dest_path = Path(output_dir) / Path(sbml_path).name
    sbml.save_converted_file(str(dest_path))
    
    constrains_file = Path(output_dir) / "constrains.tsv"
    with constrains_file.open("w") as f:
        sbml.contrains_on_kinetic_constants(f)
    
    # TODO: generate objective
if __name__ == "__main__":
    main()