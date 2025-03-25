SIS Adaptive Epidemic Simulation
================================

Description
-----------
This Python package provides two complementary SIS epidemic‑on‑adaptive-networks models:
• A stochastic Gillespie algorithm (gillespie.py)  
• A deterministic effective‑degree ODE model (effective_degree.py)  

It includes a ready‑to‑run script (`run_simulation.py`) which reproduces a test figure.

Prerequisites
-------------
• Python ≥3.8  
• networkx  
• numpy  
• scipy  
• matplotlib  

Installation
------------
Clone the repository and install in “editable” (development) mode:

    git clone <repository_url>
    cd <repository_directory>
    pip install -e .

Command‑Line Usage
------------------
Simply run:

    python run_simulation.py

This will:
1. Build an Erdős–Rényi graph (n=1000, m=1900).  
2. Execute one effective‑degree ODE simulation and five Gillespie simulations.  
3. Plot infected‑fraction vs. time and display it interactively.  

(Modify `run_simulation.py` to call `plt.savefig("simulation_output.png")` if you wish to save the figure.) :contentReference[oaicite:0]{index=0}

API Usage
---------
### Gillespie (stochastic) simulation

```python
from sis_adaptive.gillespie import Gillespie_SIS
import networkx as nx

G = nx.gnm_random_graph(1000, 1900)
t, avg_deg, inf_frac, M_vals = Gillespie_SIS(
    G, tau=0.6, gamma=1.0, alpha=0.05, omega=0.1,
    M_min=6, M_max=max(dict(G.degree()).values()), I_mid=600, k=0.05,
    initial_infecteds=list(range(100)), tmin=0, tmax=50
)
``` :contentReference[oaicite:1]{index=1}

### Effective‑degree (ODE) simulation

```python
from sis_adaptive.effective_degree import SIS_effective_degree_from_graph

times, S, I, M_vals, degree_vals, tv = SIS_effective_degree_from_graph(
    G, tau=0.6, gamma=1.0, alpha=0.05, omega=0.1,
    M_min=6, M_max=max(dict(G.degree()).values()), I_mid=600, k=0.05,
    initial_infecteds=list(range(100)), tmin=0, tmax=50, tcount=1001
)
``` :contentReference[oaicite:2]{index=2}

Outputs
-------
- **Gillespie_SIS** returns `(time, avg_degree_norm, infected_fraction_norm, adaptation_M)`  
- **SIS_effective_degree_from_graph** returns `(times, susceptible_fraction, infected_fraction, adaptation_M, degree_series, time_series)`

Project Structure
-----------------

run_simulation.py
setup.py
gillespie.py
effective_degree.py
utils.py
tests/


License:
MIT License


Author
Saiprasad v R ‹saiprasad5299@gmail.com›
