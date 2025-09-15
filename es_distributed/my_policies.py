import logging
import pickle
import numpy as np

logger = logging.getLogger(__name__)

# ------------------------
# Funzioni per flat vector
# ------------------------
class SetFromFlat:
    def __init__(self, var_list):
        self.var_list = var_list
        self.shapes = [v.shape for v in var_list]
        self.sizes = [v.size for v in var_list]
        self.total_size = sum(self.sizes)
        
    def __call__(self, theta):
        start = 0
        for i, v in enumerate(self.var_list):
            sz = self.sizes[i]
            v[:] = theta[start:start+sz].reshape(self.shapes[i])
            start += sz

class GetFlat:
    def __init__(self, var_list):
        self.var_list = var_list
        
    def __call__(self):
        return np.concatenate([v.ravel() for v in self.var_list])

# ------------------------
# Classe astratta Policy
# ------------------------
class Policy:
    def __init__(self, *args, **kwargs):
        self.args, self.kwargs = args, kwargs
        # In NumPy, _initialize deve restituire le variabili come lista di array
        self.scope = self._initialize(*args, **kwargs)
        self.all_variables = self.scope['variables']
        self.trainable_variables = self.scope['trainable_variables']

        self.num_params = sum(v.size for v in self.trainable_variables)
        self._setfromflat = SetFromFlat(self.trainable_variables)
        self._getflat = GetFlat(self.trainable_variables)

        logger.info('Trainable variables ({} parameters)'.format(self.num_params))
        for v in self.trainable_variables:
            logger.info('- shape:{} size:{}'.format(v.shape, v.size))
        logger.info('All variables')
        for v in self.all_variables:
            logger.info('- shape:{} size:{}'.format(v.shape, v.size))

    def _initialize(self, *args, **kwargs):
        """Da implementare nelle classi figlie"""
        raise NotImplementedError

    def save(self, filename):
        """Salva le variabili e gli argomenti in un file .npz"""
        np.savez(filename, **{f'v{i}': v for i, v in enumerate(self.all_variables)},
                 args_and_kwargs=pickle.dumps((self.args, self.kwargs)))

    @classmethod
    def Load(cls, filename, extra_kwargs=None):
        """Carica una policy da file .npz"""
        data = np.load(filename, allow_pickle=True)
        args, kwargs = pickle.loads(data['args_and_kwargs'].item())
        if extra_kwargs:
            kwargs.update(extra_kwargs)
        policy = cls(*args, **kwargs)
        for i, v in enumerate(policy.all_variables):
            v[:] = data[f'v{i}']
        return policy

    # ------------------------
    # Rollout
    # ------------------------
    def rollout(self, env, *, render=False, timestep_limit=None, save_obs=False, random_stream=None):
        env_timestep_limit = getattr(env.spec, 'max_episode_steps', 1000)
        timestep_limit = env_timestep_limit if timestep_limit is None else min(timestep_limit, env_timestep_limit)
        rews = []
        t = 0
        if save_obs:
            obs = []
        ob = env.reset()
        for _ in range(timestep_limit):
            ac = self.act(ob[None], random_stream=random_stream)[0]
            if save_obs:
                obs.append(ob)
            ob, rew, done, _ = env.step(ac)
            rews.append(rew)
            t += 1
            if render:
                env.render()
            if done:
                break
        rews = np.array(rews, dtype=np.float32)
        if save_obs:
            return rews, t, np.array(obs)
        return rews, t

    # ------------------------
    # Interfaccia per le classi figlie
    # ------------------------
    def act(self, ob, random_stream=None):
        raise NotImplementedError

    def set_trainable_flat(self, x):
        self._setfromflat(x)

    def get_trainable_flat(self):
        return self._getflat()

    @property
    def needs_ob_stat(self):
        raise NotImplementedError

    def set_ob_stat(self, ob_mean, ob_std):
        raise NotImplementedError