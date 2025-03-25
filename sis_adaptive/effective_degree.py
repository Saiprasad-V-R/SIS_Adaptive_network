import networkx as nx
import random
import numpy as np
from collections import defaultdict, Counter
from scipy.special import binom
from scipy.integrate import odeint
import os
import matplotlib.pyplot as plt

def _dSIS_effective_degree(X, t, original_shape, tau, gamma, alpha, omega, M_min, M_max, I_mid, k, M_values,deg,tv):
    

    ksq = original_shape[0] * original_shape[1]
    Ssi = X[:ksq]
    Isi = X[ksq:]

    Ssi.shape = original_shape
    Isi.shape = original_shape

    M = M_max - (M_max - M_min) / (1 + np.exp(-k * (Isi.sum() - I_mid)))
    M_values.append(M)
    tv.append(t)
        
    ISS = sum([sum([i * s * Ssi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    SS = sum([sum([s * Ssi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    ISI = sum([sum([i * i * Ssi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    SI = sum([sum([s * Isi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])

    PS = sum([sum([(M - (s + i)) * Ssi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    PI = sum([sum([(M - (s + i)) * Isi[s, i] for i in range(original_shape[1])]) for s in range(original_shape[0])])
    
    ML = sum([sum([(M -(s + i))*(Ssi[s, i] + Isi[s, i]) for i in range(original_shape[1])]) for s in range(original_shape[0])])
    MLM = sum([sum([((s + i))*(Ssi[s, i] + Isi[s, i]) for i in range(original_shape[1])]) for s in range(original_shape[0])])
    deg.append(MLM)  

    

    dSsi = np.zeros(original_shape)
    dIsi = np.zeros(original_shape)

    
    for s in range(original_shape[0]):
        for i in range(original_shape[1]):
            if s == 0 or i + 1 == original_shape[1]:
                Ssm1ip1 = 0
                Ism1ip1 = 0
            else:
                Ssm1ip1 = Ssi[s - 1, i + 1]
                Ism1ip1 = Isi[s - 1, i + 1]
            if i == 0 or s + 1 == original_shape[0]:
                Ssp1im1 = 0
                Isp1im1 = 0
            else:
                Ssp1im1 = Ssi[s + 1, i - 1]
                Isp1im1 = Isi[s + 1, i - 1]
            if i + 1 == original_shape[1]:
                Ssip1 = 0
                Isip1 = 0
            else: 
                Ssip1 = Ssi[s, i + 1]
                Isip1 = Isi[s, i + 1]
            if s == 0:
                Ssm1i = 0
                Ism1i = 0
            else:
                Ssm1i = Ssi[s - 1 , i]
                Ism1i = Isi[s - 1 , i]
            if i == 0:
                Ssim1 = 0
                Isim1 = 0
            else:
                Ssim1 = Ssi[s , i - 1]
                Isim1 = Isi[s , i - 1]
            if s + 1 == original_shape[0]:
                Ssp1i = 0
                Isp1i = 0
            else:
                Ssp1i = Ssi[s + 1 , i]
                Isp1i = Isi[s + 1 , i]
 
            dSsi[s, i] = -tau * i * Ssi[s, i] + gamma * Isi[s, i] \
                        + gamma * ((i + 1) * Ssm1ip1 - i * Ssi[s, i]) \
                        + tau * ISS * ((s + 1) * Ssp1im1 - s * Ssi[s, i]) / SS \
                        - omega * (s + i) * Ssi[s, i] - alpha*(M - (s + i)) * Ssi[s, i] \
                        + omega * ((i + 1)* Ssip1 + (s + 1) * Ssp1i) \
                        + alpha * (M - (s - 1 + i)) * (PS / ML) * Ssm1i \
                        + alpha * (M - (s + i - 1)) * (PI / ML) * Ssim1


            dIsi[s, i] = tau * i * Ssi[s, i] - gamma * Isi[s, i] \
                        + gamma * ((i + 1) * Ism1ip1 - i * Isi[s, i]) \
                        + tau * ISI * ((s + 1) * Isp1im1 - s * Isi[s, i]) / SI \
                        - omega * (s + i) * Isi[s, i] - alpha*(M - (s + i)) * Isi[s, i] \
                        + omega * ((i + 1)* Isip1 + (s + 1) * Isp1i) \
                        + alpha * (M - (s - 1 + i)) * (PS / ML) * Ism1i \
                        + alpha * (M - (s + i - 1)) * (PI / ML) * Isim1

    dSsi.shape = (original_shape[0] * original_shape[1])
    dIsi.shape = (original_shape[0] * original_shape[1])
    
    return np.concatenate((dSsi, dIsi), axis=0)
    
def SIS_effective_degree(Ssi0, Isi0, tau, gamma, alpha, omega, M_min, M_max, I_mid, k, num_nodes, tmin=0, tmax=50, tcount=1001):

    M_values = []  # List to store M
    deg = []
    tv = []
    times = np.linspace(tmin, tmax, tcount)
    original_shape = Ssi0.shape
    ksq = original_shape[0] * original_shape[1]
    Ssi0.shape = (1, ksq)
    Isi0.shape = (1, ksq)
    X0 = np.concatenate((Ssi0[0], Isi0[0]), axis=0)
    X = odeint(_dSIS_effective_degree, X0, times, args=(original_shape, tau, gamma,alpha,omega, M_min, M_max, I_mid, k, M_values,deg,tv))

    Ssi = X.T[0:ksq]
    Isi = X.T[ksq:]

    S = Ssi.sum(axis=0)/num_nodes
    I = Isi.sum(axis=0)/num_nodes


    return times, S, I, M_values,deg,tv

    
    
    
def SIS_effective_degree_from_graph(G, tau, gamma, alpha, omega, M_min, M_max, I_mid, k, initial_infecteds=None, 
                                    rho=None, tmin=0, tmax=100, tcount=1001):

                            
    Nk = Counter(dict(G.degree()).values())
    maxk = max(Nk.keys())
    S_si0 = np.zeros((maxk+1,maxk+1))
    I_si0 = np.zeros((maxk+1,maxk+1))
 
    if initial_infecteds is not None:
        status = _initialize_node_status_(G, initial_infecteds)
        for node in G.nodes():
            s = sum(1 for nbr in G.neighbors(node) if status[nbr]=='S')
            i = G.degree(node)-s
            if status[node] == 'S':
                S_si0[s][i] += 1
            else:
                I_si0[s][i] += 1
    else:  
        if rho is None:
            rho = 1./G.order()
        Nk = np.array([Nk[k] for k in range(maxk+1)])       
        for s in range(maxk+1):
            for i in range(maxk+1-s):
                binomial_result = binom(s+i,i)
                if binomial_result < float('Inf'):
                    S_si0[s,i] = (1-rho)*Nk[s+i] * binomial_result \
                                * (rho**i) * (1-rho)**s
                    I_si0[s,i] = rho*Nk[s+i] * binomial_result \
                                * (rho**i) * (1-rho)**s
                else:
                    S_si0[s,i] = 0
                    I_si0[s,i] = 0


    num_nodes = nx.number_of_nodes(G)

    return SIS_effective_degree(S_si0, I_si0, tau, gamma, alpha, omega, M_min, M_max, I_mid, k, num_nodes, tmin=tmin, tmax=tmax, tcount=tcount)
