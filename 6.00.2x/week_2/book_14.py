class Location(object):
    def __init__(self, x, y):
        """x and y are numbers"""
        self.x, self.y = x, y
    
    def move(self, deltaX, deltaY):
        """deltaX and deltaY are numbers"""
        return Location(self.x + deltaX, self.y + deltaY)
    
    def getX(self):
        return self.x

    def getY(self):
        return self.y
    
    def distFrom(self, other):
        ox, oy = other.x, other.y
        xDist, yDist = self.x - ox, self.y - oy
        return (xDist**2 + yDist**2)**0.5
    
    def __str__(self):
        return "<" + str(self.x) + ", " + str(self.y) + ">"
    
class Field(object):
    def __init__(self):
        self.drunks = {}
        
    def addDrunk(self, drunk, loc):
        if drunk in self.drunks:
            raise ValueError("Duplicate drunk")
        else:
            self.drunks[drunk] = loc
    
    def moveDrunk(self, drunk):
        if drunk not in self.drunks:
            raise ValueError("Drunk not in field")
        xDist, yDist = drunk.takeStep()
        currentLocation = self.drunks[drunk]
        # use move method of Location to get new location
        self.drunks[drunk] = currentLocation.move(xDist, yDist)
        
    def getLoc(self, drunk):
        if drunk not in self.drunks:
            raise ValueError("Drunk not in field")
        return self.drunks[drunk]
    
import random
class Drunk(object):
    def __init__(self, name):
        """Assumes name is str"""
        self.name = name
    
    def __str__(self):
        if self != None:
            return self.name
        return "Anonymous"
    
class UsualDrunk(Drunk):
    def takeStep(self):
        stepChoices = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        return random.choice(stepChoices)
    
def walk(f, d, numSteps):
    """Moves d numSteps times; returns the distance between the final location
    and the location at the start of the walk.

    Args:
        f (Field object): A field
        d (Drunk object): A drunk in the f field
        numSteps (int): number of steps for the drunk to take, numSteps >= 0
    """
    start = f.getLoc(d)
    for s in range(numSteps):
        f.moveDrunk(d)
    end = f.getLoc(d)
    return start.distFrom(end)

def simWalks(numSteps, numTrials, dClass):
    """Assumes numSteps an int >= 0, numTrials an int > 0,
    dClass a subclass of Drunk
    Simulates numTrials walks of numSteps steps each.
    Returns a list of the final distance for each trial"""
    Homer = dClass()