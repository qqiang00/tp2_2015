class DeBrujinGraph:
    def __init__(self, nodes: Iterable[str], k=21):
        pass # initialise la structure de données
    
    def __contains__(self, N: str) -> bool:
        pass # détermine si le graphe de Brujin contient le noeud N
    
    def __iter__(self) -> Iterable[str]:
        return self.nodes() # retourne un itérable sur les noeuds du graphe
    
    def load_factor() -> float:
        pass # calcule le facteur de charge de la table de hachage sous-jacente
    
    def add(self, N: str):
        pass # ajoute le noeud N au graphe
    
    def remove(self, N: str):
        pass # enlève le noeud N du graphe
    
    def nodes(self) -> Iterable[str]:
        pass # retourne un itérable sur les noeuds du graphe
    
    def predecessors(self, N: str) -> Iterable[str]:
        pass # retourne tous les prédécesseur du noeud N
    
    def successors(self, N: str) -> Iterable[str]:
        pass # retourne tous les successeurs du noeud N
