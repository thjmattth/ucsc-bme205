__author__ = 'tjm'

import random

class WeightedRandomizer:
    def __init__ (self, weights):
        self.__max = .0
        self.__weights = []
        for value, weight in weights.items ():
            self.__max += weight
            self.__weights.append ( (self.__max, value) )

    def random (self):
        r = random.random () * self.__max
        for ceil, value in self.__weights:
            if ceil > r: return value

# w = {'A': 1.0, 'B': 1.0, 'C': 18.0}
# # #or w = {'A': 5, 'B': 5, 'C': 90}
# # #or w = {'A': 1.0/18, 'B': 1.0/18, 'C': 1.0}
# # #or or or
# #
# wr = WeightedRandomizer (w)
#
# results = {'A': 0, 'B': 0, 'C': 0}
# for i in range (10000):
#     results [wr.random () ] += 1
#
# print ('After 10000 rounds the distribution is:')
# print (results)
