import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


planets = pd.read_csv("cumulative.csv", comment="#")
planets = planets[planets.koi_pdisposition == "CANDIDATE"]
print(planets.keys())
