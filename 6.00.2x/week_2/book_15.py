import random
def flip(numFlips):
    heads = 0
    for i in range(numFlips):
        if random.choice(("H", "T")) == "H":
            heads += 1
    return heads/numFlips

def flipSim(numFlipsPerTrial, numTrials):
    fracHeads = []
    for i in range(numTrials):
        fracHeads.append(flip(numFlipsPerTrial))
    mean = sum(fracHeads) / len(fracHeads)
    return mean

ans = flipSim(1000, 1)
print(ans)