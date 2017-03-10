import numpy as np
from utils import u_arrays as ua
import math
import matplotlib.pyplot as plt
import pandas as pd



def calc_grid_area(diameter, grid_length):

    arr = np.zeros((100,100))
    x, y = ua.draw_circle(50,50,(np.ceil(diameter / 2. /grid_length)).astype(int))

    arr[x,y]  = 1
    area_grid = np.sum(arr)* grid_length**2

    area_true = (math.pi/4) * diameter**2

    print('Diameter:', diameter)
    print("Grid area: ", area_grid)
    print("True area: ", area_true)
    print("Pixel across: ", len(np.unique(x)))

    return (area_grid, area_true)



def plot_error():
    df = pd.read_pickle('/users/global/cornkle/C_paper/wavelet/saves/pandas/3dmax_gt15000.p')

    scales = np.unique(df['scale'])

    ar = scales

    area_grid = []
    area_true = []

    for a in ar:
        grid, true = calc_grid_area(a,5)

        area_grid.append(grid)
        area_true.append(true)

    for at in np.unique(area_grid):


        pos = np.where(np.array(area_grid) == at)
        print(ar[pos[0]])

    plt.scatter(area_true, area_grid)
    plt.xlabel('true area of diameter')
    plt.ylabel('gridded area with 5km pixels')
    plt.plot(area_true, area_true, c='r')
    plt.show()

if __name__ == "__main__":
    plot_error()


