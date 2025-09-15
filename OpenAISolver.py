import math
import time
from typing import Dict, Any

from apricopt.solving.blackbox.BlackBox import BlackBox
from apricopt.solving.blackbox.BlackBoxSolver import BlackBoxSolver
import numpy as np

print("including openai code...")

from black_box_policy import BlackBoxPolicy


from es_distributed.optimizers import SGD, Optimizer, Adam

print("done!")

class OpenAISolver(BlackBoxSolver):
    
    def __init__(self):
        super().__init__()
        
    def solve(self, black_box: BlackBox, solver_params: Dict[str, Any]) -> tuple[
            Dict[str, float], float, float, float, float]:
        
        
        params: dict[str, Any] = solver_params["solver_params"]
        
        seed = params.get("seed")
        
        if seed is not None:
            np.random.seed(int(seed))
        
        N = int(params.get("population", 50))          # numero di perturbazioni per generazione
        
        episodes = int(params.get("episodes", 100))
        
        policy = BlackBoxPolicy(black_box)
        theta_init = black_box.get_np_array_initial_values()
        """lowers = black_box.get_optimization_parameters_lower_bounds_nparray()
        uppers = black_box.get_optimization_parameters_upper_bounds_nparray()
        # crea un array casuale di dimensione theta_init.size rispettando i bound per componente
        low = np.asarray(lowers, dtype=float)
        high = np.asarray(uppers, dtype=float)
        theta_init = np.array([np.random.uniform(l, h) for l, h in zip(low, high)])"""
        policy.set_trainable_flat(theta_init)
        
        
        sigma = 0.1     # deviazione standard del rumore
        stepsize = float(params.get("learning_rate", 1e-6)) # max(0.01, 0.01 * np.linalg.norm(theta_init))
        optimizer: Optimizer = SGD(policy, stepsize)
        
        
        best_f = -math.inf
        best_params = None
        best_reward = -math.inf
        
        start_time = time.time()
        
        for gen in range(episodes):
            theta = policy.get_trainable_flat()
            
            lowers = black_box.get_optimization_parameters_lower_bounds_nparray()
            uppers = black_box.get_optimization_parameters_upper_bounds_nparray()
            
            for i in range(theta.size):
                theta[i] = max(lowers[i], theta[i])
                theta[i] = min(uppers[i], theta[i])
            
            print(theta)
            
            # --- Step 1: campiona perturbazioni ---
            epsilons = np.random.randn(N, theta.size)  # perturbazioni gaussiane
            rewards = np.zeros(N)
            
            # --- Step 2: valuta reward per ogni perturbazione ---
            for i in range(N):
                policy.set_trainable_flat(theta + sigma * epsilons[i])
                rew, _ = policy.rollout()  # rollout ritorna np.array([reward])
                if -rew[0] >= black_box.get_objective_upper_bound():
                    print("[FATAL ERROR] simulator crashed")
                    exit(1)
                if rew[0] > best_f:
                    print(f"New minimum!!: {-rew[0]}")
                    best_f = rew[0]
                    best_params = theta + sigma * epsilons[i]
                rewards[i] = rew[0]  # ES usa reward scalare
                print(rewards[i])
            if best_f >= -params.get("global_minimum", -math.inf):
                break

            # --- Step 3: calcola gradiente stimato ---
            
            #result = np.zeros(theta.size)
            #for (rew, eps) in zip(rewards, epsilons):
            #    result += rew*eps
                
            #globalg = result*(1/N*sigma)
            
            # A = (rewards - np.mean(rewards)) / (np.std(rewards) + 1e-8)
            # ATTENZIONE: ho tolto il learning rate da qua
            globalg = (np.dot(epsilons.T, rewards) / (N * sigma))
            
            policy.set_trainable_flat(theta + stepsize*globalg)
            # --- Step 4: aggiorna theta usando SGD ---
            # ratio = optimizer.update(globalg)

            # --- Logging ---
            best_local = np.mean(rewards)
            if best_reward > best_local:
                print("[WARNING] This generation performed worse than the previous one")
            best_reward = best_local
            print(f"Gen {gen:3d} | mean reward: {best_reward:.6f} | Update ratio: {stepsize:.6f}")

        # theta finale
        theta_final = policy.get_trainable_flat()
        print("Miglior theta trovato suggeriuto:", theta_final)
        print(f"Migliore theta trovato: {best_params}, f: {best_f}")
        
        param_dict: Dict[str, float] = dict()
        params_ids = black_box.get_optimization_parameters_ids()
        for i in range(len(params_ids)):
            param_dict[params_ids[i]] = float(best_params[i])
            
        print(f"time: {time.time()-start_time}")
        
        return param_dict, float(best_f), 0, None, None