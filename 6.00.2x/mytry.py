# def song_playlist(songs, max_size):
#     """
#     songs: list of tuples, ('song_name', song_len, song_size)
#     max_size: float, maximum size of total songs that you can fit

#     Start with the song first in the 'songs' list, then pick the next 
#     song to be the one with the lowest file size not already picked, repeat

#     Returns: a list of a subset of songs fitting in 'max_size' in the order 
#              in which they were chosen.
#     """
#     # maximize sum of song_len, with constraint sum of song_size <= max_size
#     if len(songs) == 0:
#         return []
#     if songs[0][2] > max_size:
#         return []
#     else:
#         # Make a copy of the songs list and initialize remaining size variable
#         songs_arr = songs[:]
#         rem_size = max_size
#         # Add the first song to chosen list, update remaining size, remove that song from the array
#         chosen_arr = [songs_arr[0][0]]
#         rem_size -= songs_arr[0][2]
#         songs_arr.pop(0)
#         # Fill the secondary songs
#         while rem_size >= 0:
#             songsize_list = [song[2] for song in songs_arr]
#             mask = [song[2] == min(songsize_list) for song in songs_arr]
#             chosen = [song for song, is_included in zip(songs_arr, mask) if is_included]
#             if chosen == []:
#                 return chosen_arr
#             rem_size -= chosen[0][2]
#             if rem_size >= 0:
#                 chosen_arr.append(chosen[0][0])
#                 songs_arr.remove(chosen[0])
        
#         return chosen_arr
    
# # songs = [('Roar',4.4, 4.0),('Sail',3.5, 7.7),('Timber', 5.1, 6.9),('Wannabe',2.7, 1.2)]
# # max_size = 12.2
# # songs = [('Roar',4.4, 4.0),('Sail',3.5, 7.7),('Timber', 5.1, 6.9),('Wannabe',2.7, 1.2)]
# # max_size = 11

# # print(song_playlist(songs, max_size))

# print(song_playlist([('a', 4.4, 4.0), ('b', 3.5, 7.7), ('c', 5.1, 6.9), ('d', 2.7, 1.2)], 20))


# def max_contig_sum(L):
#     """ L, a list of integers, at least one positive
#     Returns the maximum sum of a contiguous subsequence in L """
#     sums = []
#     for idx, element in enumerate(L):
#         sums += [sum(L[idx:i+1]) for i in range(idx, len(L))]
#     return max(sums)

# # print(max_contig_sum([3, 4, -1, 5, -4]))
# print(max_contig_sum([3, 4, -1, 5, -4]))
# print(max_contig_sum([3, 4, -8, 15, -1, 2]))


def solveit(test):
    """ test, a function that takes an int parameter and returns a Boolean
        Assumes there exists an int, x, such that test(x) is True
        Returns an int, x, with the smallest absolute value such that test(x) is True 
        In case of ties, return any one of them. 
    """
    x = 0
    if test(x):
        return x
    
    x_p = 0
    x_n = 0
    while not (test(x_p) or test(x_n)):
        x_p += 1
        x_n -= 1
        
    if test(x_p):
        return x_p
    else:
        return x_n

#### This test case prints 49 ####
def f(x):
    return (x+15)**0.5 + x**0.5 == 15
print(solveit(f))

#### This test case prints 0 ####
def f(x):
    return x == 0
print(solveit(f))