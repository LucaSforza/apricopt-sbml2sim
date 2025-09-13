from typing import Set
from pandas import DataFrame

from apricopt.model.Parameter import Parameter

class MyParameter(Parameter):
    
        def __init__(self, param_id: str, name: str, parameter_scale="lin", lower_bound: float = float("-inf"),
                 upper_bound: float = float("inf"), nominal_value: float = float("-inf"),
                 distribution: str = 'uniform', mu: float = 0, sigma: float = 0,
                 granularity: float = 0, optimize: bool = True):
            super().__init__(param_id, name, lower_bound, upper_bound, nominal_value, distribution, mu, sigma, granularity, optimize)
            self.parameter_scale = parameter_scale

def my_get_parameter_space(parameter_df: DataFrame) -> Set[Parameter]:
    params: Set[Parameter] = set()
    for param_id, data in parameter_df.iterrows():
        params.add(MyParameter(str(param_id), data.parameterName,
                             data.parameterScale, data.lowerBound, data.upperBound, data.nominalValue,
                             data.distribution, data.mu, data.sigma, data.granularity, data.estimate))
    return params