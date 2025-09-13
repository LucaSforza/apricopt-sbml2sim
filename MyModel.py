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

import copy
from typing import Dict, Set, List
from apricopt.model.FastObservable import FastObservable
from apricopt.model.LogObservable import LogObservable
from apricopt.model.ModelInstance import ModelInstance
from apricopt.model.Observable import Observable, compute_expressions_values
from apricopt.model.Model import Model
from apricopt.model.Parameter import Parameter
from MyParam import MyParameter


class MyModel(Model):
    """
    Objects of this class represents models of systems that can be simulated and that may be subject to
    simulation-based optimization.

    A Model object may have an objective function and one or more constraints.
    Constraints are subdivided into regular constraints and fast constraints.
    A regular constraint is a constraint whose value can only be computed via simulation of the modes. A fast constraint,
    on the other side, does not require a simulation to be computed.
    A Model object contains the definition of its parameter space, as a set of Parameter objects.
    An object of this class has a reference to a SimulationEngine object that actually simulates it and to a
    ModelInstance object that represents the actual object that is used to simulate the system.
    In some sense, the ModelObject provides a wrapping over the ModelInstance object, which may use other frameworks to
    perform simulations and other tasks.
    Finally, the object contains further information, i.e. the integration tolerances,
    the time step, the name of the model outputs that are actually observedand the name of the file that contains the
    model (if any).
    """

    def __init__(self, sim_engine, model_filename: str, abs_tol: float,
                 rel_tol: float, time_step: float, observed_outputs: List[str] = None,
                 observed_outputs_observable: List[Observable] = None, **exec_info):
        self.sim_engine = sim_engine
        self.model_filename: str = model_filename
        if sim_engine:
            self.instance: ModelInstance = sim_engine.load_model(model_filename, **exec_info)
            self.instance.set_simulation_configuration(abs_tol, rel_tol, time_step)
            self.initial_values = self.build_initial_values()
        self.absolute_tolerance = abs_tol
        self.relative_tolerance = rel_tol
        self.time_step = time_step
        self.parameters: Dict[str, Parameter] = dict()
        self.objective: Observable = None
        self.constraints: Set[Observable] = set()
        self.fast_constraints: Set[FastObservable] = set()
        self.log_observables: Set[LogObservable] = set()
        self.cached_ids = None
        self.observed_outputs: List[str] = list() if observed_outputs is None else observed_outputs
        self.observed_outputs_observable: List[Observable] = set() if observed_outputs_observable is None \
                                                                                    else observed_outputs_observable

    def reinit(self):
        """
        Re-initializes the state of the object of type ModelInstance in the <pre>instance</pre> field.
        :return:
        """
        self.instance.__delete__(self)
        self.instance = self.sim_engine.load_model(self.model_filename)
        self.instance.set_simulation_configuration(self.absolute_tolerance, self.relative_tolerance, self.time_step)
        self.initial_values = self.build_initial_values()

    def _new_model_instance(self) -> None:
        """
        Initializes the <pre>instancce</pre> field with a ModelInstance object built starting from the model_filename.
        :return: None
        """
        self.instance = self.sim_engine.load_model(self.model_filename)
        self.instance.set_simulation_configuration(self.absolute_tolerance, self.relative_tolerance, self.time_step)
        self.initial_values = self.build_initial_values()

    def copy(self):
        """
        Performs a copy of self. This new object is a shallow copy of self, but the model instance of the copy is a
        new model instance of the same system.
        :return: None
        """
        copied_model = copy.copy(self)
        copied_model._new_model_instance()
        return copied_model

    def deep_copy(self):
        #TODO: returns a completely different Model object
        pass

    def set_parameter_space(self, params: Set[Parameter]) -> None:
        """
        Updates the parameter space using the Parameter objects in the given <pre>params</pre> set.
        :param params: a set of Parameter object. Each Parameter will be added to the parameter space of self.
        :return: None
        """
        for param in params:
            self.parameters[param.id] = param

    def set_params(self, params_values: Dict[str, float]) -> None:
        for param_id, value in params_values.items():
            param_obj = self.parameters[param_id]
            if param_obj.distribution == "uniform":
                if value < param_obj.lower_bound or \
                        value > param_obj.upper_bound:
                    raise ValueError(f"The value {value} for parameter {param_obj.name} is out of bounds.")
            if isinstance(param_obj, MyParameter):
                if param_obj.parameter_scale == "log":
                    self.instance.set_parameter(param_id, 10**params_values[param_id])
                elif param_obj.parameter_scale == "lin":
                    self.instance.set_parameter(param_id, params_values[param_id])
                else:
                    raise ValueError()
            else:
                # assumi lineare
                self.instance.set_parameter(param_id, params_values[param_id])
        self.initial_values = self.build_initial_values()

    def set_parameters_nominal_values(self, params_values: Dict[str, float]) -> None:
        for param_id, value in params_values.items():
            param_obj = self.parameters[param_id]
            if param_obj.distribution == "uniform":
                if value < param_obj.lower_bound or \
                        value > param_obj.upper_bound:
                    raise ValueError(f"The value {value} for parameter {param_obj.name} is out of bounds.")
            param_obj.nominal_value = value

    def build_initial_values(self) -> Dict[str, float]:
        return dict(self.instance.get_parameters_initial_values(),
                    **self.instance.get_compartment_initial_volumes())

    def set_fixed_params(self, params_values: Dict[str, float]) -> None:
        # for param_id, value in params_values.items(): # commented because inefficient: one single set of all parameters is much quicker than many set of one parameter
        # self.instance.set_parameter(param_id, value)
        # self.initial_values[param_id] = value         WRONG: does not consider inter-parameter dependencies
        self.instance.set_parameters(params_values)
        self.initial_values = self.build_initial_values()

    def evaluate_fast_constraints(self, params_values: Dict[str, float]) -> Dict[str, float]:
        result: Dict[str, float] = dict()
        full_assignment = self.complete_assignment(params_values)
        for fast_constraint in self.fast_constraints:
            result[fast_constraint.id] = fast_constraint.evaluate(full_assignment)
        return result

    def build_zero_sim_output(self) -> Dict[str, float]:
        sim_output: Dict[str, float] = dict()

        if self.objective:
            sim_output[self.objective.id] = self.objective.upper_bound

        for constraint in self.constraints:
            sim_output[constraint.id] = constraint.upper_bound
         
        ''''#: ToCheck    
        for feasibility_constraint in self.feasibility_constraints:
            sim_output[feasibility_constraint.id] = feasibility_constraint.upper_bound'''
        return sim_output

    def get_observables_ids(self) -> Set[str]:
        if not self.cached_ids:
            self.cached_ids: Set[str] = set()
            if self.objective:
                self.cached_ids.add(self.objective.id)

            for constraint in self.constraints:
                self.cached_ids.add(constraint.id)

        return self.cached_ids

    def evaluate_constraints(self, trajectory: Dict[str, List[float]]) -> Dict[str, float]:
        full_trajectory: Dict[str, List[float]] = self.complete_trajectory(trajectory)
        result: Dict[str, float] = dict()
        if self.objective:
            result[self.objective.id] = self.objective.evaluate(full_trajectory)
        for constraint in self.constraints:
            result[constraint.id] = constraint.evaluate(full_trajectory)
        
        '''# ToCheck
        for feasibility_constraint in self.feasibility_constraints:
            result[feasibility_constraint.id] = feasibility_constraint.evaluate(full_trajectory)
        '''
        return result

    def evaluate_log_observables(self, trajectory: Dict[str, List[float]]) -> Dict[str, str]:
        full_trajectory: Dict[str, List[float]] = self.complete_trajectory(trajectory)
        result: Dict[str, str] = dict()
        for log_observable in self.log_observables:
            result[log_observable.id] = log_observable.evaluate(full_trajectory)
        return result

    def complete_trajectory(self, trajectory: Dict[str, List[float]]) -> Dict[str, List[float]]:
        length = len(trajectory['time'])
        full_trajectory: Dict[str, List[float]] = dict(trajectory)
        for value_id, value in self.initial_values.items():
            if value_id not in full_trajectory:
                full_trajectory[value_id] = [value] * length

        return full_trajectory

    def complete_assignment(self, assignment: Dict[str, float]) -> Dict[str, float]:
        full_assignment: Dict[str, float] = dict(assignment)
        for value_id, value in self.initial_values.items():
            if value_id not in full_assignment:
                full_assignment[value_id] = value
        return full_assignment

    def compute_first_satisfaction_index(self, obs: Observable, trajectory: Dict[str, List[float]]) -> int:
        # assume that function is Identity, so only one expression
        if len(obs.expressions) > 1:
            raise ValueError("An initialization constraint cannot have more than one expression in its formula")

        full_trajectory: Dict[str, List[float]] = self.complete_trajectory(trajectory)
        constraint_values = compute_expressions_values(full_trajectory, obs.expressions, obs.function)[0]
        trajectory_length: int = len(full_trajectory['time'])

        for t in range(trajectory_length):
            if constraint_values[t] <= 0:
                return t
        return trajectory_length
