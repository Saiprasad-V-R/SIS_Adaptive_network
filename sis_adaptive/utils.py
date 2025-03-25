def initialize_node_status(G, initial_infecteds):
    """
    Return a dict mapping each node in G to 'I' if infected, else 'S'.
    """
    return {n: ('I' if n in initial_infecteds else 'S') for n in G.nodes()}

