import matplotlib.pyplot as plt
import sys
import re
import atexit


def main():
    log_file_path = sys.argv[1]
    attempts = []
    best_losses = []
    with open(log_file_path, "r") as f:
        for line in f.readlines():
            if line.startswith("Attempt:"):
                m = re.search(r'Attempt:\s*(\d+),\s*loss\s*([^,]+),\s*best loss\s*([^\s,]+)', line)
                if not m:
                    print(f"[FATAL ERROR] regular expression is invalid: {m}")
                    exit(1)
                attempt = int(m.group(1))
                _loss = float(m.group(2))
                best_loss = float(m.group(3))

                attempts.append(attempt)
                best_losses.append(best_loss)
    plt.figure()
    plt.plot(attempts, best_losses, linestyle='-')
    plt.yscale('log')
    plt.xlabel('Attempt')
    plt.ylabel('Best loss')
    plt.title('Ottimizzazione')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig("plot.png")

if __name__ == "__main__":
    main()