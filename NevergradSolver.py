import math
import time
from typing import Dict, Any

from apricopt.solving.blackbox.BlackBox import BlackBox
from apricopt.solving.blackbox.BlackBoxSolver import BlackBoxSolver

import numpy as np

class NevergradSolver(BlackBoxSolver):
    import nevergrad as ng
    import numpy as np

    class Permutation(ng.p.Array):
        def __init__(self, choices: list[Any], **kwargs):
            self._choices = choices
            super().__init__(shape=(len(choices),), lower=0.0, upper=1.0, **kwargs)

        @property
        def value(self):
            # decodifica in permutazione di 0..n-1
            indices = np.argsort(super().value)
            return [self._choices[i] for i in indices]


    def __init__(self, verbose = False):
        super().__init__()
        self.verbose = verbose
    
    def solve(self, black_box: BlackBox, solver_params: Dict[str, Any]) -> tuple[
            Dict[str, float], float, float, float, float]:
        import nevergrad as ng
        
        parametrization: Dict[str, ng.p.Parameter] = {}
        
        for param in black_box.get_optimization_parameters_ids():
            lower = black_box.get_optimization_parameter_lower_bound(param)        
            upper = black_box.get_optimization_parameter_upper_bound(param)
            init = black_box.get_optimization_parameter_initial_value(param)
            gran = black_box.get_optimization_parameter_granularity(param)
            if gran:
                # TODO: add init
                parametrization[param] = ng.p.TransitionChoice(np.arange(lower, upper+gran, gran))
            else:
                parametrization[param] = ng.p.Scalar(lower=lower, upper=upper, init=init)        
        parametrization: ng.p.Dict = ng.p.Dict(**parametrization)
        seed = solver_params.get("seed")
        if seed is not None:
            parametrization.random_state.seed(int(seed))
        constrains_path: str = solver_params.get("constrains")
        if constrains_path:
            pairs: list[tuple[str,str]] = []
            with open(constrains_path, "r") as f:
                for line in f.readlines():
                    parameters = line.split()
                    controller = parameters[0]
                    controlled = parameters[1]
                    pairs.append((controller, controlled))
            def constraint(x: ng.p.Dict) -> float:
                loss = 0.0
                for controller, controlled in pairs:
                    value = x.value[controlled] - x.value[controller]
                    if value > 0:
                        loss += value
                return -loss
            parametrization.register_cheap_constraint(constraint)
        solver_parameters = solver_params["solver_params"]
        budget: int = solver_parameters.get("budget", 3000)
        optimizer_name: str = solver_parameters.get("optimizer", "NGOpt")
        
        optimizer_cls = getattr(ng.optimizers, optimizer_name, ng.optimizers.NGOpt)
        if optimizer_name == "BayesianOptimization":
            opt = ng.optimizers._BO(parametrization=parametrization, budget=budget,init_budget=20, utility_kind="ei")
        else:
            opt = optimizer_cls(parametrization=parametrization, budget=budget, num_workers=1)
        print("Nevergrad optimzier chosen: "+ getattr(optimizer_cls, "__name__", type(optimizer_cls).__name__))
        best_value = math.inf
        best_parameters: Dict[str, float] = None
        start_time = time.time()
        for i in range(budget):
            params = opt.ask()
            losses = black_box.evaluate(params.value)
            loss = 0.0
            for (_,value) in losses.items():
                loss += value
            if loss < best_value:
                best_value = loss
                best_parameters = params.value
            if best_value <= 0:
                break
            opt.tell(params, loss)
            if self.verbose:
                print(f"Attempt: {i}, loss {loss}, best loss {best_value}")
        print(f"[INFO] ended search in {time.time() - start_time}")
        return best_parameters, best_value, 0, None, None