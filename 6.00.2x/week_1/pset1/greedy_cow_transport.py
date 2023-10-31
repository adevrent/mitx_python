def greedy_cow_transport(cows, wlim):
    """Greedy algorithm for selection of transporting cows with spaceship

    Args:
        cows (dict): a dict of cows where keys are cow names (str) and values are cow weights (int)
        wlim (int): weight limit of the spaceship
    """
    cows_sorted_list = sorted(list(cows.items()), key = lambda x: x[1], reverse = True)
    total_transport = []
    while len(cows_sorted_list) != 0:
        wcount = 0
        to_transport = []
        for name, weight in cows_sorted_list[:]:
            if wcount + weight <= wlim:
                to_transport.append(name)
                wcount += weight
        total_transport.append(to_transport)
        cows_sorted_list = [tup for tup in cows_sorted_list[:] if tup[0] not in to_transport]
    return total_transport

mytrans = greedy_cow_transport({'Milkshake': 75, 'Polaris': 20, 'Louis': 45, 'Miss Bella': 15, 'Patches': 60, 'Muscles': 65, 'Horns': 50, 'Clover': 5, 'MooMoo': 85, 'Lotus': 10}, 100)
print(mytrans)