import networkx as nx
import random
import numpy as np
from collections import defaultdict, Counter
from scipy.special import binom
from scipy.integrate import odeint
import os
import matplotlib.pyplot as plt

def Gillespie_SIS(G, tau, gamma, alpha, omega, M_min, M_max, I_mid, k, initial_infecteds=None, tmin=0, tmax=100):
    """
    We would like to express our sincere gratitude to the authors of 'Mathematics of Epidemics on Networks' 
    for their invaluable contributions to the field of epidemiology. 
    Their work has provided the foundational framework for our research and codebase.
    Simulate an SIS model with dynamic network adaptation using the Gillespie algorithm.

    Parameters:
    - G (networkx.Graph): The initial graph structure representing the network.
    - tau (float): Transmission rate of infection.
    - gamma (float): Recovery rate of infection.
    - alpha (float): Link creation rate.
    - omega (float): Link removal rate.
    - M_min (float): Minimum adaptation parameter.
    - M_max (float): Maximum adaptation parameter.
    - I_mid (float): Infection midpoint for M calculation.
    - k (float): Steepness of the M function.
    - initial_infecteds (list[int]): List of initially infected node indices.
    - tmin (float): Start time for simulation.
    - tmax (float): End time for simulation.

    Returns:
    - t_values (np.ndarray): Array of time points in the simulation.
    - k_values (np.ndarray): Array of normalized average degrees over time.
    - c_values (np.ndarray): Array of normalized infection counts over time.
    - M_values (np.ndarray): Array of adaptation parameter M values over time.
    """

    # Initialize rates, statuses, and tracking lists
    tau = float(tau)
    gamma = float(gamma)

    I = [len(initial_infecteds)]
    S = [G.order() - I[0]]
    times = [tmin]
    t_values = []
    k_values = []
    M_values = []
    c_values = []  # Store values of c
    
    infecteds = set(initial_infecteds)
    IS_links = set()


    status = defaultdict(lambda: 'S')
    for node in initial_infecteds:
        status[node] = 'I'
        for nbr in G.neighbors(node):
            if status[nbr] == 'S':
                IS_links.add((node, nbr))

    t = tmin

    # Calculate initial M
    M = M_max 
    t_values.append(t)
    c_values.append(I[0])
    M_values.append(M)
    
    # Find the number of nodes in the graph
    num_nodes = nx.number_of_nodes(G)
    average_degree = sum(dict(G.degree()).values())/2  
    nmx = ((num_nodes *M) - average_degree*2)/2
    k_values.append(average_degree)

    total_recovery_rate = gamma * len(infecteds)
    total_transmission_rate = tau * len(IS_links)
    total_link_removal_rate = omega * average_degree 
    total_link_creation_rate = alpha * nmx
    
    total_rate = total_recovery_rate + total_transmission_rate + total_link_removal_rate + total_link_creation_rate
    delay = random.expovariate(total_rate)
    t = t+delay

    edges_list = list(G.edges())

    while  t < tmax:
        c = 0
        for node in G.nodes():
            if status[node] =='I':
                c+=1
                for nbr in G.neighbors(node):
                    if status[nbr] == 'S':
                        IS_links.add((node, nbr))

        c_values.append(c)  # Store the value of c
        
        
        
        if random.random() < total_recovery_rate / total_rate:  # recover
            recovering_node = random.choice(list(infecteds))
            infecteds.remove(recovering_node)
            status[recovering_node] = 'S'  # Update the status of the recovering node
            
            for nbr in G.neighbors(recovering_node):
                if nbr == recovering_node:  # move past self edges
                    continue
                elif status[nbr] == 'S' and G.has_edge(recovering_node, nbr):
                    IS_links.discard((recovering_node, nbr))
                else:
                    IS_links.add((nbr, recovering_node))
            
            times.append(t)
            S.append(S[-1] + 1)
            I.append(I[-1] - 1)


        if random.random() < total_transmission_rate / total_rate:  # transmit
                transmitter, recipient = random.choice(list(IS_links))
                if (transmitter, recipient) in G.edges() or (recipient, transmitter) in G.edges():
                    status[recipient] = 'I'
                    infecteds.add(recipient)
            
                    for nbr in G.neighbors(recipient):
                        if status[nbr] == 'S':
                            IS_links.add((recipient, nbr))
                        elif nbr != recipient and ((nbr, recipient)) in IS_links:
                            IS_links.remove((nbr, recipient))

                    times.append(t)
                    S.append(S[-1] - 1)
                    I.append(I[-1] + 1)
        
        
        if random.random() < total_link_removal_rate  /  total_rate: #link removal
            if G.number_of_edges() > 1:
                random_node = random.choice(list(G.nodes()))
                degree = G.degree(random_node)
                edges_list = list(G.edges())
                link_to_remove = random.choice(edges_list)
                G.remove_edge(*link_to_remove)
                if link_to_remove in IS_links:
                    IS_links.remove(link_to_remove)

        
        if random.random() < total_link_creation_rate / total_rate:  # link creation
            nodes = random.sample(list(G.nodes()), 2)
            node1, node2 = nodes
            edge = (node1, node2)
            G.add_edge(*edge)
            node1_status = status[node1]
            node2_status = status[node2]
            if node1_status == 'I' and node2_status == 'S':
                IS_links.add(edge)
            elif node1_status == 'S' and node2_status == 'I':
                IS_links.add((node2, node1))  # Ensure the order is (infected, susceptible)

        M = M_max - ((M_max - M_min) / (1 + np.exp(-k * (c - I_mid))))
        
        M_values.append(M)
        average_degree = sum(dict(G.degree()).values())/2  
        nmx = ((num_nodes * M) - average_degree*2)/2
        k_values.append(average_degree)
        
        total_recovery_rate = gamma * len(infecteds)
        total_transmission_rate = tau * len(IS_links)
        total_link_removal_rate = omega * average_degree
        total_link_creation_rate = alpha * nmx
        total_rate = total_recovery_rate + total_transmission_rate + total_link_removal_rate + total_link_creation_rate
        
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
        t_values.append(t)
        
    return np.array(t_values),np.array(k_values)/num_nodes, np.array(c_values)/num_nodes,np.array(M_values)  # Return c_values along with other outputs





