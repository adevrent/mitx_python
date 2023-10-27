class Node(object):
    def __init__(self, name) -> None:
        """Assumes name is a string"""
        self.name = name
    def getName(self):
        return self.name
    def __str__(self) -> str:
        return self.name

class Edge(object):
    def __init__(self, src, dest) -> None:
        """Assumes src and dest are nodes"""
        self.src = src
        self.dest = dest
    def getSource(self):
        return self.src
    def getDestination(self):
        return self.dest
    def __str__(self) -> str:
        return self.src.getName() + "->" + self.dest.getName()
    
class WeightedEdge(Edge):
    def __init__(self, src, dest, weight):
        # Your code here
        self.src = src
        self.dest = dest
        self.weight = weight
    def getWeight(self):
        # Your code here
        return self.weight
    def __str__(self):
        # Your code here
        return str(self.src) + "->" + str(self.dest) + " " + "(" +  str(self.weight) + ")"