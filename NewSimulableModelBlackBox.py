"""
This file is part of Apricopt.

Apricopt is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Apricopt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Apricopt.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2020-2021 Marco Esposito, Leonardo Picchiami.
"""
from typing import List, Dict

import numpy as np

from apricopt.model.Model import Model
from apricopt.simulation.SimulationEngine import SimulationEngine
from apricopt.simulation.roadrunner.RoadRunnerEngine import RoadRunnerEngine
from apricopt.solving.blackbox.BlackBox import BlackBox
import json

import sbml2sim.s2s as s2s

class NewSimulableModelBlackBox(BlackBox):

    def __init__(self, sim_engine: SimulationEngine, model: Model, horizon: float, loss_function: s2s.ErrorHandler):
        if sim_engine.model_instance_class() != model.instance.__class__:
            raise TypeError("The provided SimulationEngine object must be able to simulate the model.")
        self.sim_engine = sim_engine
        self.model = model
        self.horizon = horizon
        self.param_ids = None
        self.eb_constraint_ids = None
        self.pb_constraint_ids = None
        self.loss_function = loss_function
        self._granularity_is_required=False

    def evaluate(self, parameters: Dict[str, float], check_input=True) -> Dict[str, float]:
        # ME: Force parameters within bounds
        for p_id, p_value in parameters.items():
            if p_value < self.get_optimization_parameter_lower_bound(p_id):
                parameters[p_id] = self.get_optimization_parameter_lower_bound(p_id)
            elif p_value > self.get_optimization_parameter_upper_bound(p_id):
                parameters[p_id] = self.get_optimization_parameter_upper_bound(p_id)
        # Set parameters in the model
        self.model.set_params(parameters)
        # Simulate the model using its simulation engine
        result: Dict[str, float] = dict()
        
        try:
            trajectory: dict[str, list[float]] = self.sim_engine.simulate_trajectory(self.model, self.horizon)
            # Ensure all values are JSON-serializable (convert numpy arrays and numpy numbers)
            """serializable = {}
            for k, v in trajectory.items():
                if isinstance(v, np.ndarray):
                    serializable[k] = v.tolist()
                else:
                    serializable[k] = [float(x) for x in v]
            print(json.dumps(serializable, indent=4, ensure_ascii=False))"""
        except:
            result[self.get_objective_id()] = self.get_objective_upper_bound()
            return result
        # Extract the results from the trajectory produced by the simulator
        

        cost = self.loss_function.from_trajectory(trajectory, 0.5) # TODO: parametrize

        if cost > self.get_objective_upper_bound():
            cost = self.get_objective_upper_bound() 

        result[self.get_objective_id()] = cost

        # Non ci sono vincoli, commento
        # for ebc_id in self.get_extreme_barrier_constraints_ids():
        #     result[ebc_id] = trajectory[ebc_id][-1]
        # for pbc_id in self.get_progressive_barrier_constraints_ids():
        #     result[pbc_id] = trajectory[pbc_id][-1]
        return result

    def evaluate_np_array(self, params: np.array, check_input=False) -> Dict[str, float]:
        param_dict: Dict[str, float] = dict()
        params_ids = self.get_optimization_parameters_ids()
        for i in range(len(params_ids)):
            param_dict[params_ids[i]] = params[i]
        return self.evaluate(param_dict)

    def evaluate_objective_np_array(self, params: np.array, check_input=False) -> float:
        param_dict: Dict[str, float] = dict()
        params_ids = self.get_optimization_parameters_ids()
        for i in range(len(params_ids)):
            param_dict[params_ids[i]] = params[i]
        return self.evaluate_objective(param_dict)

    def is_input_valid(self, parameters: Dict[str, float]) -> bool:
        return True

    def get_optimization_parameters_number(self) -> int:
        return len(self.model.parameters)

    def get_optimization_parameters_ids(self) -> List[str]:
        if not self.param_ids:
            p_ids = [param.id for param in self.model.parameters.values()]
            p_ids.sort()
            self.param_ids = p_ids
        return self.param_ids

    def get_optimization_parameter_lower_bound(self, param_id) -> float:
        return self.model.parameters[param_id].lower_bound

    def get_optimization_parameter_upper_bound(self, param_id) -> float:
        return self.model.parameters[param_id].upper_bound

    def get_optimization_parameters_lower_bounds_nparray(self) -> np.array:
        lb = []
        for param_id in self.get_optimization_parameters_ids():
            lb += [self.model.parameters[param_id].lower_bound]
        return np.array(lb)

    def get_optimization_parameters_upper_bounds_nparray(self) -> np.array:
        ub = []
        for param_id in self.get_optimization_parameters_ids():
            ub += [self.model.parameters[param_id].upper_bound]
        return np.array(ub)

    def get_optimization_parameter_initial_value(self, param_id) -> float:
        return self.model.parameters[param_id].nominal_value

    # TODO implement caching
    def get_optimization_parameters_initial_values(self) -> np.array:
        res: Dict[str, float] = dict()
        for param_id in self.get_optimization_parameters_ids():
            res[param_id] = self.model.parameters[param_id].nominal_value
        return res

    def get_optimization_parameters_initial_values_np_array(self) -> np.array:
        nv = []
        for param_id in self.get_optimization_parameters_ids():
            nv += [self.model.parameters[param_id].nominal_value]
        return np.array(nv)

    def set_optimization_parameters_initial_values(self, param_values: Dict[str, float]) -> None:
        self.model.set_parameters_nominal_values(param_values)

    def optimization_parameters_initial_values_are_empty(self) -> bool:
        return False

    def get_optimization_parameter_granularity(self, param_id) -> float:
        return self.model.parameters[param_id].granularity

    def get_extreme_barrier_constraints_number(self) -> int:
        return len(self.model.fast_constraints)

    def get_progressive_barrier_constraints_number(self) -> int:
        return len(self.model.constraints)

    def get_extreme_barrier_constraints_ids(self) -> List[str]:
        if not self.eb_constraint_ids:
            c_ids = [constraint.id for constraint in self.model.fast_constraints]
            c_ids.sort()
            self.eb_constraint_ids = c_ids
        return self.eb_constraint_ids

    def get_progressive_barrier_constraints_ids(self) -> List[str]:
        if not self.pb_constraint_ids:
            c_ids = [constraint.id for constraint in self.model.constraints]
            c_ids.sort()
            self.pb_constraint_ids = c_ids
        return self.pb_constraint_ids

    def get_objective_id(self) -> str:
        return "cost"

    def get_objective_upper_bound(self) -> float:
        return 1e60

    @staticmethod
    def get_raisable_exception_type():
        return TypeError

    def finalize(self) -> None:
        pass

    def granularity_is_required(self) -> bool:
        return self._granularity_is_required

    def set_granularity_is_required(self, val: bool) -> None:
        self._granularity_is_required = val
