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
