
# coding: utf-8

# In[34]:


import numpy as np
from itertools import combinations
from scipy.sparse import dok_matrix
from operator import add

class SimplicialComplex:
    def __init__(self, simplices=[]):
        self.import_simplices(simplices=simplices)
        print('1')
 
    def import_simplices(self, simplices=[]):
        self.simplices = map(lambda simplex: tuple(sorted(simplex)), simplices)
        #self.face_set = self.faces()
        self.face_set = faces(self)


# In[35]:


def faces(self):
    faceset = set()
    for simplex in self.simplices:
        numnodes = len(simplex)
        for r in range(numnodes, 0, -1):
            for face in combinations(simplex, r):
                faceset.add(face)
    return faceset


# In[36]:


import networkx as nx
from scipy.spatial import distance
from itertools import product
 
class VietorisRipsComplex(SimplicialComplex):
    def __init__(self, points, epsilon, labels=None, distfcn=distance.euclidean):
        self.pts = points
        self.labels = range(len(self.pts)) if labels==None or len(labels)!=len(self.pts) else labels
        self.epsilon = epsilon
        self.distfcn = distfcn
        self.network = self.construct_network(self.pts, self.labels, self.epsilon, self.distfcn)
        self.import_simplices(map(tuple, list(nx.find_cliques(self.network))))
 
    def construct_network(self, points, labels, epsilon, distfcn):
        g = nx.Graph()
        g.add_nodes_from(labels)
        zips = zip(points, labels)
        for pair in product(zips, zips):
            if pair[0][1]!=pair[1][1]:
                dist = distfcn(pair[0][0], pair[1][0])
                if dist < epsilon:
                    g.add_edge(pair[0][1], pair[1][1])
        return g


# In[37]:


def n_faces(self, n):
    return filter(lambda face: len(face)==n+1, self.face_set)


# In[38]:


print('2')
sampl = np.random.uniform(low=-6.0, high=6.0, size=(50,3))
#sc = SimplicialComplex([('a', 'b', 'c', 'd')])
sc = SimplicialComplex(sampl)
print(sc)
faces(sc)
#sg=VietorisRipsComplex(sc)
sg2=AlphaComplex(sc,'1')


# In[ ]:


from scipy.spatial import Delaunay, distance
from operator import or_
from functools import partial
 
def facesiter(simplex):
    for i in range(len(simplex)):
        yield simplex[:i]+simplex[(i+1):]
 
def flattening_simplex(simplices):
    for simplex in simplices:
        for point in simplex:
            yield point

def get_allpoints(simplices):
    return set(flattening_simplex(simplices))
 
def contain_detachededges(simplex, distdict, epsilon):
    if len(simplex)==2:
        return (distdict[simplex[0], simplex[1]] > 2*epsilon)
    else:
        return reduce(or_, map(partial(contain_detachededges, distdict=distdict, epsilon=epsilon), facesiter(simplex)))
 
class AlphaComplex(SimplicialComplex):
    def __init__(self, points, epsilon, labels=None, distfcn=distance.euclidean):
        self.pts = points
        self.labels = range(len(self.pts)) if labels==None or len(labels)!=len(self.pts) else labels
        self.epsilon = epsilon
        self.distfcn = distfcn
        self.import_simplices(self.construct_simplices(self.pts, self.labels, self.epsilon, self.distfcn))
 
    def calculate_distmatrix(self, points, labels, distfcn):
        distdict = {}
        for i in range(len(labels)):
            for j in range(len(labels)):
                distdict[(labels[i], labels[j])] = distfcn(points[i], points[j])
        return distdict
 
    def construct_simplices(self, points, labels, epsilon, distfcn):
        delaunay = Delaunay(points)
        delaunay_simplices = map(tuple, delaunay.simplices)
        distdict = self.calculate_distmatrix(points, labels, distfcn)
 
        simplices = []
        for simplex in delaunay_simplices:
            faces = list(facesiter(simplex))
            detached = map(partial(contain_detachededges, distdict=distdict, epsilon=epsilon), faces)
            if reduce(or_, detached):
                if len(simplex)>2:
                    for face, notkeep in zip(faces, detached):
                        if not notkeep:
                            simplices.append(face)
            else:
                simplices.append(simplex)
        simplices = map(lambda simplex: tuple(sorted(simplex)), simplices)
        simplices = list(set(simplices))
 
        allpts = get_allpoints(simplices)
        for point in (set(labels)-allpts):
            simplices += [(point,)]
 
        return simplices

