#!/usr/bin/python
import pickle

parameters = pickle.load(open("parameters.pickle", "rb"))
print(str(parameters))
