from typing import Dict, Any

from apricopt.solving.blackbox.BlackBox import BlackBox
from apricopt.solving.blackbox.BlackBoxSolver import BlackBoxSolver
import numpy as np

print("including openai code...")

from black_box_policy import BlackBoxPolicy


from es_distributed.optimizers import SGD, Optimizer

print("done!")

class OpenAISolver(BlackBoxSolver):
    
    def __init__(self):
        super().__init__()
        
    def solve(self, black_box: BlackBox, solver_params: Dict[str, Any]) -> tuple[
            Dict[str, float], float, float, float, float]:
        policy = BlackBoxPolicy(black_box)
        theta_init = black_box.get_np_array_initial_values()
        policy.set_trainable_flat(theta_init)
        N = 50          # numero di perturbazioni per generazione
        sigma = 0.1     # deviazione standard del rumore
        stepsize = 0.01 * np.linalg.norm(theta_init)
        optimizer: Optimizer = SGD(policy, stepsize)
        
        params: dict[str, Any] = solver_params["solver_params"]
        
        budget = int(params.get("budget", 100))
        optimizer.update()
        
        for gen in range(budget):
            theta = policy.get_trainable_flat()
            
            # --- Step 1: campiona perturbazioni ---
            epsilons = np.random.randn(N, theta.size)  # perturbazioni gaussiane
            rewards = np.zeros(N)
            
            # --- Step 2: valuta reward per ogni perturbazione ---
            for i in range(N):
                policy.set_trainable_flat(theta + sigma * epsilons[i])
                rew, _ = policy.rollout()  # rollout ritorna np.array([reward])
                rewards[i] = rew[0]  # ES usa reward scalare

            # --- Step 3: calcola gradiente stimato ---
            A = (rewards - np.mean(rewards)) / (np.std(rewards) + 1e-8)
            globalg = np.dot(epsilons.T, A) / (N * sigma)

            # --- Step 4: aggiorna theta usando SGD ---
            ratio = optimizer.update(globalg)

            # --- Logging ---
            best_reward = np.max(rewards)
            print(f"Gen {gen:3d} | Best reward: {best_reward:.6f} | Update ratio: {ratio:.6f}")

        # theta finale
        theta_final = policy.get_trainable_flat()
        print("Miglior theta trovato:", theta_final)