__all__ = ["centergrid"]

import numpy as np
# INPUT: "array": Array of N values asumed evenly spaced in lin or log space
#        "spacing": string either "linear" or "log" specifying spacing of array
# OUTPUT: 
#        "array_new": Array of length N+1, containing linear or log midpoints of values in "array",
#                     extrapolated one gridpoint in each end of the midpoint array. An array of
#                     of midpoints would have length N-1, so exptrapolating 1 point in each end
#                     gives length N+1.

def centergrid(array,spacing):
    if spacing == 'linear':
        array_new = (array[:-1]+array[1:])/2.
        array_new = np.insert(array_new,0,array[0]-(array_new[0]-array[0]))
        array_new = np.append(array_new,array[-1] + (array[-1]-array_new[-1]))
    elif spacing == 'log':
        array_new = np.sqrt(array[:-1]*array[1:]) 
        array_new = np.insert(array_new, 0, array_new[0]  / (array_new[1] /array_new[0]))
        array_new = np.append(array_new,    array_new[-1] * (array_new[-1]/array_new[-2]))
    return array_new

