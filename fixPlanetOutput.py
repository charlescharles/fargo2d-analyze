import numpy as np


if __name__ == '__main__':
    sec = np.loadtxt('bigplanet0.dat')
    date = sec[:, -2]
    _, ixs = np.unique(date, return_index=True)
    newSec = sec[ixs]
    np.save('sec', newSec)