# intoexo

intoexo (INTeriors Of EXOplanets) is meant to simulate the interiors of planets using simple spherical models
that incorporate a small amount of chemistry.

## planet.py
Where the implementation of the model is using scipy.quad to integrate (my first attempt)

# planetdfq.py
Uses DFQs and the Euler method to integrate. I think this is how its supposed to be done.

## analysis.ipynb
A notebook meant to interface with planet.py and use it to create plots and analyze the model.
