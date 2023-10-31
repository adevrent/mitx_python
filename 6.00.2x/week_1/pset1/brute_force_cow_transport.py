from ps1_partition import *
def brute_force_cow_transport(cows,limit=10):
    """
    Finds the allocation of cows that minimizes the number of spaceship trips
    via brute force.  The brute force algorithm should follow the following method:

    1. Enumerate all possible ways that the cows can be divided into separate trips
    2. Select the allocation that minimizes the number of trips without making any trip
        that does not obey the weight limitation
            
    Does not mutate the given dictionary of cows.

    Parameters:
    cows - a dictionary of name (string), weight (int) pairs
    limit - weight limit of the spaceship (an int)
    
    Returns:
    A list of lists, with each inner list containing the names of cows
    transported on a particular trip and the overall list containing all the
    trips
    """
    total_partitions = [prt for prt in get_partitions(cows.items())]
    candidate_journeys = []
    for journey in total_partitions:
        flag = True
        for trip in journey:
            if sum([tup[1] for tup in trip]) > limit:
                flag = False
                break
        if flag:
            candidate_journeys.append(journey)
    
    return sorted(candidate_journeys, key = lambda x: len(x))[0]