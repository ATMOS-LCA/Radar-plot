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

from datetime import datetime

#A partir da lista datalist, cria uma sub-lista apenas com os dados desejados
def typeFilter(datalist):
    typeData = input(f'Informe o tipo de dado: ')

    filteredData = []
    for data in datalist:
        indice_dot = data.find('.')
        if typeData == data[16:indice_dot]:
            filteredData.append(data)

    return filteredData
    #print(filteredData)

def timeFilter(datalist):
    data_inicial_str = input(f'Informe uma data inicial (dd/mm/yyyy): ')
    data_final_str = input(f'Informe uma data final (dd/mm/yyyy): ')
    
    #Converte as datas de entrada para o formato datetime a partir do formato dd/mm/yyyy
    data_inicial = datetime.strptime(data_inicial_str, "%d/%m/%Y")
    data_final = datetime.strptime(data_final_str, "%d/%m/%Y")

    filteredData = []
    for data in datalist:
        #Converte a data do arquivo para o formato datetime a partir do formato yyyymmdd
        data_arquivo = datetime.strptime(data[0:8], "%Y%m%d")
        if data_inicial <= data_arquivo <= data_final:
            filteredData.append(data)

    return filteredData
    #print(filteredData)

#Acessa a pasta data e cria uma lista com os arquivos que estão lá
datalist = os.listdir('./data')
arquivos = timeFilter(typeFilter(datalist))
for arquivo in arquivos:
    fpath = './data/' + arquivo
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

    gif_path = './gifs/'+arquivo+'.gif'
    imageio.mimsave(gif_path, frames, fps=1, loop=0)
