import numpy as np
data = np.load('lasagne4bio/subcellular_localization/data/train.npz')
for varName in data:
    print data[varName][:]

    
