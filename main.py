import json

import apricopt.IO.data_input
# from apricopt.simulation.COPASI.COPASIEngine import COPASIEngine
from apricopt.simulation.SimulationEngine import SimulationEngine
from apricopt.simulation.roadrunner.RoadRunnerEngine import RoadRunnerEngine
from apricopt.solving.blackbox.BlackBoxSolver import BlackBoxSolver
from apricopt.solving.blackbox.pyswarms.PySwarmsSolver import PySwarmsSolver
from apricopt.solving.blackbox.scipy.SciPySolver import SciPySolver

from apricopt.solving.blackbox.NOMAD.NOMADSolver import NOMADSolver
from NevergradSolver import NevergradSolver
from NewSimulableModelBlackBox import NewSimulableModelBlackBox
from apricopt.IO.data_input import *
import petab
import random
print("[INFO] inclued apricopt")
from MyModel import MyModel

import sbml2sim.s2s as s2s

from MyParam import my_get_parameter_space

'''
Input da linea di comando:
1. id del paziente
    da questo si "calcola" il path del suo config.yaml
    si calcola anche la directory dove scrivere l'output
2. Numero di run di NOMAD

'''
import sys

path_to_dir_example = sys.argv[1]

config_filename = f"{path_to_dir_example}/config.yaml"
N = 10 # TODO: cosa fa N?
'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading input configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

# Read the YAML file
data: dict = apricopt.IO.data_input.parse_config_file(config_filename)

# Instantiate the simulation engine
sim_engine: SimulationEngine
#if data['simulator'].lower() == "copasi":
#    sim_engine = COPASIEngine()
#el
if data['simulator'].lower() == "roadrunner":
    sim_engine = RoadRunnerEngine()
else:
    raise ValueError()

# Instantiate the solver
solver: BlackBoxSolver
if data['solver'].lower() == "nomad":
    solver = NOMADSolver()
elif data['solver'].lower() == "pyswarms":
    solver = PySwarmsSolver()
elif data['solver'].lower() == "scipy":
    solver = SciPySolver()
elif data['solver'].lower() == "nevergrad":
    solver = NevergradSolver(verbose=True)
elif data['solver'].lower() == "openai":
    from OpenAISolver import OpenAISolver
    print("[INFO] inclued openai")
    solver = OpenAISolver()
else:
    raise ValueError()

# Read other info from the YAML file
files = data['files']
abs_tol = data['absolute_tolerance']
rel_tol = data['relative_tolerance']
time_step = data['time_step']
horizon = float(data['horizon'])
random_seed = int(data['random_seed'])
loss_functions: dict[str, dict[str, int]] = data['loss_functions']

# If present, read the observed outputs from the json file
if 'observed_outputs' in files:
    f = open(files['observed_outputs'])
    observed_outputs_text = f.read()
    f.close()
    observed_outputs = json.loads(observed_outputs_text)
else:
    observed_outputs = None

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating and configuring the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

# Load the SBML model and instantiate it.
model = MyModel(sim_engine, f"{path_to_dir_example}/{files['model']}", abs_tol, rel_tol, time_step, observed_outputs=observed_outputs)

# Create the PEtab problem. This helps in reading the files and validates the information contained in them
problem = petab.problem.Problem.from_files(sbml_file=f"{path_to_dir_example}/{files['model']}",
                                           parameter_file=f"{path_to_dir_example}/{files['parameters']}",
                                           observable_files=f"{path_to_dir_example}/{files['objective']}")

# Create the parameter space from the PEtab parameter table
parameter_space: Set[Parameter] = my_get_parameter_space(problem.parameter_df)
# Set the parameter space for the optimization problem
model.set_parameter_space(parameter_space)

# Superfluo, ma lo lasciamo
# Create the objective observable from the PEtab observable table
# objective: Observable = get_objective(problem.observable_df)
# Set the Observable for the objective function
# model.objective = objective

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating and configuring the black-box
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

error_sum = s2s.error_sum_create()
for (loss_f_name, settings) in loss_functions.items():
    scale = float(settings["scale"])
    function_name = f"{loss_f_name}_error_create"
    handler = getattr(s2s, function_name, None)
    if handler is None or not callable(handler):
        raise AttributeError(f"Function '{function_name}' not found in sbml2sim.s2s")
    loss_f = None
    if loss_f_name == "fitting":
        file = f"{path_to_dir_example}/{files['ordinamento']}"
        loss_f = s2s.fitting_error_create()
        with open(file) as f:
            for line in f.readlines():
                less = line.split()[0]
                more = line.split()[1]
                s2s.fitting_error_add_constrain(loss_f, less, more)
    elif loss_f_name == "value":
        file = f"{path_to_dir_example}/{files['valori']}"
        loss_f = s2s.value_error_create()
        with open(file) as f:
            for line in f.readlines():
                id = line.split()[0]
                value = float(line.split()[1])
                s2s.value_error_add_constrain(loss_f, id, value)
    else:
        loss_f = handler()
    s2s.error_sum_add_handler(error_sum, loss_f, scale)

# Instantiate the black box
black_box = NewSimulableModelBlackBox(sim_engine, model, horizon, error_sum)

'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
# Get the solver hyperparameters from the YAML file
solver_parameters = data['solver_parameters']

print(f"\n\n======= Optimization ========\n")
constrains =files.get("constrains")
path_to_constrains = None 
if constrains:
    path_to_constrains = f"{path_to_dir_example}/{constrains}"

# Create the hyper-parameter map
solver_params = {"solver_params": solver_parameters}
if path_to_constrains and data['solver'].lower() == "nevergrad":
    solver_parameters["constrains"] = path_to_constrains
print(solver_params)

# Run the optimization
result = solver.solve(black_box, solver_params)

# Extract info from the result
optimal_params, best_objective, h_return, n_of_bb_evaluations, n_of_iterations = result
print(f"\n\n======= Results ========\n")
print(result)
