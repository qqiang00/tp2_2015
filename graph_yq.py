import random
from collections.abc import MutableMapping
from random import randrange

def random_seq(bases = "ATCG", k = 21):
    return "".join(random.choices(bases, k = k))

'''consider a DNA sequence with fixed length 21:
there are at most: pow(4, 21) = 4^21 = 4398046511104 possible sequences.
'''

def seq_hash_code(seq: str) -> int:
    """convert a dna sequence (as a key) to a hash code
    params
        seq: dna sequence with bases from: 'A','T','C','G'
    returns
        code: hash code for key k, int 
    """
    bases = "ATCG"
    codes = ["1"] # highest bit, if len(seq) is quarantee to fixed 21
                  # we may not use highest bit 1.
    for base in seq:
        # use two bits represent each base: A:00, T:01, C:10, G:11
        codes.append(bin(bases.index(base))[2:].zfill(2)) # binary code

    code = int("".join(codes), 2) # converto binary code to decimal
    return code


# MyProbHashMap class using MAD (multiply-add-and-divide, MAD) as compression 
# function and linear probing collision-handling scheme
class MyProbeHashMap(MutableMapping):
    class _Item:
        __slots__ = '_key', '_value'
        def __init__( self, k, v = None ):
            self._key = k
            self._value = v

        def __eq__( self, other ):
            return self._key == other._key

        def __ne__( self, other ):
            return not( self == other )

        def __lt__( self, other ):
            return self._key < other._key

        def __ge__( self, other ):
            return self._key >= other._key

        def __str__( self ):
            return "<" + str( self._key ) + "," + str( self._value ) + ">"

        def key( self ):
            return self._key

        def value( self ):
            return self._value    


    _AVAIL = object()
    _MIN_CAP = 11               # minimal capacity of the Map
    
    def __init__( self, cap = 11, p = 109345121, probing = 'linear' ):
        cap = max(MyProbeHashMap._MIN_CAP, cap)  # minimal capacity of Map
        # is this prime p big enough for k = 21?
        self._T = cap * [None]               # bucket array with cap entries
        self._n = 0                          # number of elements in array
        self._prime = p                      # prime number for MAD compression
        # a: scale
        self._scale = 1 + randrange( p - 1 ) # scale(also called: a) in [1, p-1]
        # seems we don't need this loop as p is prime
        trouve = False
        while not trouve:
            self._scale = 1 + randrange( p - 1 )
            if not ( self._scale % p ) == 0:
                trouve = True
        
        self._shift = randrange( p )         # shift(also called: b) in [0, p-1]
        self._mask = cap                     # N
        self._q = 7                          # used for double hash
        self._prob_func = self._probing_linear
        if probing == "quadratic":
            self._prob_func = self._probing_quadratic
        elif probing == "double_hash":
            self._prob_func = self._probing_double_hash
        elif probing == "random":
            self._prob_func = self._probing_random
            
        self._collisions = 0


    # hash function using _hash_code and MAD compression function
    def _hash_function( self, k ):
        tmp = seq_hash_code(k) * self._scale + self._shift 
        tmp = tmp  % self._prime % self._mask
        return tmp        
        
    def get_load_factor(self):
        return self._n / len(self._T)
        
    # length (number of elements) in the table(T)
    def __len__( self ):
        return self._n

    def _is_available( self, j ):
        return self._T[j] is None or self._T[j] is MyProbeHashMap._AVAIL
    

    # three different probing methods to solve collision
    def _probing_linear(self, i, *args):
        return i
    
    def _probing_quadratic(self, i, *args):
        """probing quadratic"""
        return i*i
    
    def _probing_double_hash(self, i, k, *args):
        """double hashing
        f(i) = i * h'(k)
        h'(k) = q - (k mod q)
        """
        h_prime = self._q - (seq_hash_code(k) % self._q)
        return i * h_prime

    def _probing_random(self, i, k, *args):
        # this probing doesn't work well.
        random.seed(1)
        return i * random.randint(1, 7)
    
    def get_collisions(self):
        return self._collisions
    
    def reset_collisions(self):
        self._collisions = 0
        
    def _find_slot(self, j, k):
        """compared with _find_slot, this method support three different
        probing methods, and we can count the collision numbers
        """
        first_avail = None
        i = 0
        while True:
            new_j = (j + self._prob_func(i, k)) % len(self._T)
            if self._is_available(new_j):# May be None or _AVAIL
                if first_avail is None:
                    first_avail = new_j
                if self._T[new_j] is None:
                    return (False, first_avail)
            elif k == self._T[new_j]._key:
                return (True, new_j)
            # collision:
            i += 1
            self._collisions += 1
                
            
    def _bucket_getitem( self, j, k ):
        found, s = self._find_slot( j, k )
        if not found:
            return False, None
        return True, self._T[s]._value


    def _bucket_setitem( self, j, k, v ):
        found, s = self._find_slot( j, k )
        if not found:
            self._T[s] = self._Item( k, v )
            self._n += 1
        else:
            self._T[s]._value = v


    def _bucket_delitem( self, j, k ):
        found, s = self._find_slot( j, k )
        if not found:
            return False, None
        value = self._T[s]._value
        self._T[s] = MyProbeHashMap._AVAIL
        return True, value   # what if a value itself is "False"?

    def __iter__( self ):
        for j in range( len( self._T ) ):
            if not self._is_available( j ):
                yield self._T[j]._key

    def __getitem__( self, k ):
        j = self._hash_function( k )
        success, value = self._bucket_getitem( j, k )
        if not success: # do not use value to check the condition,
                        # what if value itself is False
            raise KeyError
        else:
            return value

    def __setitem__( self, k, v ):
        j = self._hash_function( k )
        self._bucket_setitem( j, k, v )
        if self._n > len( self._T ) // 2:
            self._resize( 2 * len( self._T ) - 1 )

    def __delitem__( self, k ):
        j = self._hash_function( k )
        succes, _ = self._bucket_delitem( j, k )
        if succes: # do not use value to check the condition,
                   # what if value itself is False
            self._n -= 1
            # need resize? minimal size = 
            if self._n < len(self._T) // 4:
                new_cap = max(MyProbeHashMap._MIN_CAP, (len(self._T)+1)/2)
                self._resize(new_cap)
        else:
            raise KeyError

    def _resize( self, c ):
        old = list( self.items() )
        self._T = c * [None]
        self._n = 0
        self._mask = c
        for (k,v) in old:
            self[k] = v

    def is_empty( self ):
        return len( self ) == 0

    def __str__( self ):
        if self.is_empty():
            return "{}"
        pp = "{"
        for item in self.items():
            pp += str( item )
        pp += "}"
        pp += " size = "
        pp += str( len( self ) )
        return pp

    def get( self, k, d = None ):
        try:
            tmp = self[k]
            return tmp
        except KeyError:
            return d

    def setdefault( self, k, d = None ):
        try:
            tmp = self[k]
            return tmp
        except:
            self[k] = d
            return d


class DeBrujinGraph:
    _bases = "ACTG"
    class Edge:
        __slots__ = '_origin', '_destin', '_element'
        #constructeur avec élément
        def __init__( self, u, v, x ):
            self._origin = u
            self._destin = v
            self._element = x

        #retourne les Sommets origine et destination dans un tuple
        def endpoints( self ):
            return( self._origin, self._destin )

        #retourne le Sommet opposé à v
        def opposite( self, v ):
            return self._destin if v is self._origin else self._origin

        #accès à l'élément
        def element( self ):
            return self._element

        #prettyprint d'une Arête
        def __str__( self ):
            return str( self._element )

        #la référence au tuple (origin, destination) comme hashcode
        def __hash__( self ):
            #hachage du tuple (origin, destination)
            return seq_hash_code(self._origin + self._destin)

    def __init__(self, nodes, k = 21):
        self._outgoing = MyProbeHashMap()
        self._incoming = MyProbeHashMap()
        self._k = k
        for node in nodes: # insert all node into graph
            self._insert_node(node)
        # build edges
        for node in self.nodes():
            successors = self._all_possible_successors(node)
            for successor in successors:
                if successor in self.nodes():
                    self._insert_edge(node, successor, successor[-1])
                    
    
    def __str__(self):
        s = "G( nodes{ "
        for v in self.nodes():
            s += str( v ) + " "
        s += "}, edge{ "
        for e in self._edges():
            s += str( e ) + " "
        s += "} )"
        return s

    def _all_possible_successors(self, N: str) -> [str]:
        return [N[1:]+base for base in DeBrujinGraph._bases]

    def _all_possible_predecessors(self, N: str) -> [str]:
        return [base + N[:-1] for base in DeBrujinGraph._bases]
        

    def _edges(self):
        result = set()
        for secondary_map in self._outgoing.values():
            result.update(secondary_map.values())
        return result    
  

    def _validate(self, x, edge = False):
        if x is None:
            raise ValueError("x is None")
        if not isinstance(x, str):
            raise TypeError("x should be an instance of str, now is:{}".format(type(x)))
        if edge and (len(x) != 1 or x not in DeBrujinGraph._bases):
            raise ValueError("x should be one of character in ACTG")
        if not edge and (len(x) != self._k):
            raise ValueError("x should have {} characters".format(self._k))

    def _insert_node(self, x = None):
        self._validate(x)
        if x in self.nodes():
            return  x # x already in graph
        self._outgoing[x] = MyProbeHashMap()
        self._incoming[x] = MyProbeHashMap()
        return x

    def _insert_edge(self, u, v, x = None):
        self._validate(u)
        self._validate(v)
        self._validate(x, True)
        edge = self.Edge(u, v, x)
        self._outgoing[u][v] = edge
        self._incoming[v][u] = edge


    def __contains__(self, N: str) -> bool:
        # détermine si le graphe de Brujin contient le noeud N
        return not self._outgoing.get(N, None)

    def load_factor(self) -> float:
        # calcule le facteur de charge de la table de hachage sous-jacente
        return self._outgoing.get_load_factor()
        
    def add(self, N: str):
        # ajoute le noeud N au graphe
        self._validate(N)
        if N in self.nodes():
            print("{} already exits in graph".format(N))
            return 
        self._insert_node(N)
        successors = self._all_possible_successors(N)
        for successor in successors:
            if successor in self.nodes():
                self._insert_edge(N, successor, successor[-1])
        predecessors = self._all_possible_predecessors(N)
        for predecessor in predecessors:
            if predecessor in self.nodes():
                self._insert_edge(predecessor, N, N[-1])
    
    def remove(self, N: str):
        self._validate(N)
        if not N in self.nodes():
            print("deleting... {} not in graph".format(N))
            return
        successors = self._all_possible_successors(N)
        for successor in successors:
            if successor in self.nodes():
                del self._outgoing[N][successor]
        predecessors = self._all_possible_predecessors(N)
        for predecessor in predecessors:
            if predecessor in self.nodes():
                del self._incoming[N][predecessor]
        del self._outgoing[N]
        del self._incoming[N]
        #del self._incoming[N]
        pass # enlève le noeud N du graphe
    
    def nodes(self):
        # retourne un itérable sur les noeuds du graphe
        return self._outgoing.keys()
    
    def predecessors(self, N: str):
        # retourne tous les prédécesseur du noeud N
        self._validate(N)
        all_possible = self._all_possible_predecessors(N)
        result = []
        for predecessor in all_possible:
            if predecessor in self.nodes():
                result.append(predecessor)
        return result
        
    
    def successors(self, N: str):
        # retourne tous les successeurs du noeud N
        self._validate(N)
        all_possible = self._all_possible_successors(N)
        result = []
        for successor in all_possible:
            if successor in self.nodes():
                result.append(successor)
        return result


if __name__ == "__main__":

    seq = 'ACTGAGTC'
    k = 2
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]

    print(kmers)
    print(seq)

    graph = DeBrujinGraph(kmers, k=2)

    print(graph)

    print('************')
    for node in graph.nodes():
        print(node, 'successor = ',graph.successors(node), "predecessors =",graph.predecessors(node) )


    print('load factor =' , graph.load_factor())

    ##########################################
    print('\r\n','### add GA ###')
    graph.add("GA")
    for node in graph.nodes():
        print(node, 'successor = ',graph.successors(node), "predecessors =",graph.predecessors(node) )

    print('\r\n','### remove GA ###')
    graph.remove("GA")
    for node in graph.nodes():
        print(node, 'successor = ',graph.successors(node), "predecessors =",graph.predecessors(node) )

