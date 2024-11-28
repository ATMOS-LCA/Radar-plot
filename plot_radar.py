import xarray, h5py
import xradar as xd
import wradlib as wrl
import numpy as np
import pylab as pl
import os
import imageio

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyart

#Acessa a pasta data e cria uma lista com os arquivos que estão lá
datalist = os.listdir('./data')

#A partir da lista datalist, cria uma sub-lista apenas com os dados desejados
typeData = input(f'Informe o tipo de dado: ')

especificData = []
for data in datalist:
    if typeData in data:
        especificData.append(data)

#Tranforma da data de dd/mm/yyyy para yyyymmdd
day, month, year = map(str, input(f'Informe a data que deseja visualizar: ').split('/'))
Day = year + month + day

fpath = './data/' + datalist[0]
frames = []
data = wrl.io.hdf.read_generic_hdf5(os.path.join(fpath))
for i in range(1, 14):
    fig = plt.figure(figsize=(10,10))
    da = wrl.georef.create_xarray_dataarray(data["dataset"+str(i)+"/data1/data"]["data"]).wrl.georef.georeference()
    im = da.wrl.vis.plot(fig=fig, crs="cg")

    aux = "dataset"+str(i)
    frame_path = f'./frames/{aux}.png'
    plt.savefig(frame_path, dpi=100)
    frames.append(imageio.imread(frame_path))
    plt.close()

gif_path = './gifs/dBZ.gif'
imageio.mimsave(gif_path, frames, fps=1)