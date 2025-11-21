import networkx as nx
from sis_adaptive.gillespie import Gillespie_SIS
from sis_adaptive.effective_degree import SIS_effective_degree_from_graph
import matplotlib.pyplot as plt
from collections import defaultdict, Counter


initial_size = 100
n = 1000  # Number of nodes
m = 1900  # Number of edges
tmax = 50
rho = initial_size / n
tau = 0.6  # Transmission rate
gamma = 1  # Recovery rate
initial_infected = range(initial_size)
omega = 0.1
alpha = 0.05
iterations = 5
M_min = 6
I_mid = 200
k_value = 0.05  
G = nx.gnm_random_graph(n, m)

# Calculate max degree for M_max
Nk = Counter(dict(G.degree()).values())
maxk = max(Nk.keys())
M_max = maxk

# Run SIS effective degree simulation
t_sis, S_sis, I_sis, M_values, deg, tv = SIS_effective_degree_from_graph(
    G, tau, gamma, alpha, omega, M_min, M_max, I_mid, k_value, rho=0.1, tmax=tmax
)

# Plot the SIS effective degree results
plt.plot(t_sis, I_sis, linewidth=4, label='Effective Degree Model')

# Run Gillespie simulations and plot each iteration
iterations = 5
for counter in range(iterations):
    G1 = G.copy()
    t, S, I, Mm = Gillespie_SIS(
        G1, tau, gamma, alpha, omega, M_min, M_max, I_mid, k_value, initial_infected, tmin=0, tmax=tmax
    )
    if counter == 0:
        plt.plot(t, I, color='k', alpha=0.09, label='Simulation')
    else:
        plt.plot(t, I, color='k', alpha=0.09)

# Final plot adjustments
plt.xlabel('Time')
plt.ylabel('Infected Population')
plt.title('SIS Model with Effective Degree vs Gillespie Simulation')
plt.legend()
plt.show()
