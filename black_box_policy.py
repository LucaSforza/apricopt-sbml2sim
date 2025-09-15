import numpy as np

print("including policies...")

from es_distributed.my_policies import Policy

from apricopt.solving.blackbox.BlackBox import BlackBox

class BlackBoxPolicy(Policy):
    def __init__(self, black_box: BlackBox):
        self.black_box = black_box
        self.num_params = black_box.get_optimization_parameters_number()
        # dummy per compatibilità con ES
        self._setfromflat = lambda x: setattr(self, 'theta', x)
        self._getflat = lambda: getattr(self, 'theta', np.zeros(self.num_params))

    def rollout(self, env=None, **kwargs):
        # ES chiama rollout e si aspetta (rews, t)
        theta = self.get_trainable_flat()
        reward_dict = self.black_box.evaluate_np_array(theta)
        # Convertiamo il dizionario in reward scalare
        reward = -reward_dict[self.black_box.get_objective_id()]  # ES massimizza
        return np.array([reward]), 1  # rews, timesteps

    def act(self, ob, random_stream=None):
        raise NotImplementedError  # non serve

    @property
    def needs_ob_stat(self):
        return False

    def set_ob_stat(self, ob_mean, ob_std):
        pass