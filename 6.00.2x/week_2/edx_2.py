class Location(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def getX(self):
        return self.x
    
    def getY(self):
        return self.y
    
    def move(self, deltaX, deltaY):
        return Location(self.x + deltaX, self.y + deltaY)
    
    def distFrom(self, other):
        diffX = self.x - other.x
        diffY = self.y - other.y
        dist = (diffX**2 + diffY**2) ** 0.5
        return dist
    
    def __str__(self):
        return "<" + str(self.x) + ", " + str(self.y)


class Field(object):
    def __init__(self):
        self.drunks = {}
        
    def addDrunk(self, drunk, loc):
        if drunk in self.drunks:
            raise ValueError("Duplicate drunk")
        self.drunks[drunk] = loc
        
    def getLoc(self, drunk):
        if drunk not in self.drunks:
            raise ValueError("Drunk not in field")
        return self.drunks[drunk]
    
    def moveDrunk(self, drunk):
        if drunk not in self.drunks:
            raise ValueError("Drunk not in field")
        xDist, yDist = drunk.takeStep()
        currentLocation = self.drunks[drunk]
        newLocation = currentLocation.move(xDist, yDist)
        self.drunks[drunk] = newLocation


class Drunk(object):
    def __init__(self, name):
        self.name = name
        
    def __str__(self):
        return "This drunk is named " + str(self.name)
    
    
import random
class UsualDrunk(Drunk):
    def takeStep(self):
        stepChoices = [(0.0, 1.0), (0.0, -1.0), (1.0, 0.0), (-1.0, 0.0)]
        return random.choice(stepChoices)


class ColdDrunk(Drunk):
    def takeStep(self):
        stepChoices = [(0.0, 0.9), (0.0, -1.1), (1.0, 0.0), (-1.0, 0.0)]
        return random.choice(stepChoices)
    
    
def walk(f, d, numSteps):
    """Simulate a walk consisting of numSteps number of steps.
    The action includes a single Drunk, which is in Field f.

    Args:
        f (Field class instance): a Field object
        d (Drunk class instance): a Drunk object which is in f
        numSteps (int): Number of steps to take in the walk. numSteps >= 0
    """
    startingLoc = f.getLoc(d)
    for step in range(numSteps):
        f.moveDrunk(d)
    endLoc = f.getLoc(d)
    return endLoc.distFrom(startingLoc)


def simWalks(numSteps, numTrials, dClass):
    """Simulates numTrials number of walks of a Drunk subclass (UsualDrunk or ColdDrunk),
    saves the distance difference of each trial to a list named distances.
    Returns distances.

    Args:
        numSteps (int): # of steps for each walk
        numTrials (int): # of trials total in the simulation
        dClass (Drunk subclass): A subclass of Drunk
    """
    distances = []
    origin = Location(0, 0)
    for t in range(numTrials):
        f = Field()
        d = dClass("d")
        f.addDrunk(d, origin)
        distances.append(walk(f, d, numSteps))
    return distances


def drunkTest(walkLengths, numTrials, dClass):
    """Assumes walkLengths a sequence of ints >= 0
         numTrials an int > 0, dClass a subclass of Drunk
       For each number of steps in walkLengths, runs simWalks with
         numTrials walks and prints results"""
    for numSteps in walkLengths:
        distances = simWalks(numSteps, numTrials, dClass)
        print(dClass.__name__, 'random walk of', numSteps, 'steps')
        print(' Mean =', round(sum(distances)/len(distances), 4))
        print(' Max =', max(distances), 'Min =', min(distances))
        

random.seed(0)
drunkTest((10, 100, 1000, 10000), 100, UsualDrunk)