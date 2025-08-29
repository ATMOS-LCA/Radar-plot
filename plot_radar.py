#!apt-get update
#!apt-get install -y libgdal-dev libhdf5-dev libnetcdf-dev

# Instalacao de bibliotecas essenciais para radar e dados
#!pip install xarray h5py numpy scipy matplotlib pylab imageio
#!pip install xradar
#!pip install wradlib
#!pip install pyart
#!pip install plotly
#!pip install -U kaleido # O Kaleido eh necessario para exportar graficos como imagens estaticas

# Instalacao de bibliotecas para geoprocessamento e mapeamento
#!pip install cartopy
#!pip install shapely # Geralmente instalada com cartopy, mas bom garantir
#!pip install geopy
#!pip install geographiclib

import xarray as xr
#from open_radar_data import DATASETS
import xradar as xd

import h5py
import wradlib as wrl
import numpy as np
import pylab as pl
import os
import imageio.v2 as imageio

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import pyart
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from geopy.distance import geodesic
from geographiclib.geodesic import Geodesic

from datetime import datetime

import plotly.graph_objects as go
import pandas as pd

#A partir da lista fornecida, cria uma sub-lista apenas com os dados desejados
def typeFilter(datalist):
    typeData = input(f'Informe o tipo de dado: ')

    filteredData = []
    for data in datalist:
      indice_dot = data.find('.')
      if typeData == data[16:indice_dot]:
        filteredData.append(data)

    return filteredData
    #print(filteredData)

#A partir da lista fornecida, cria uma sub-lista apenas com o período de tempo especificado
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

#Codigo para gerar o plot a partir dos arquivos .vol.h5
def plotRadar(datalist):
  data_vol = []
  #th_vals = []
  for fname in datalist:
      fpath = '/content/drive/MyDrive/data/' + fname #Alterar o caminho da pasta, caso nao seja executado no colab
      data = wrl.io.hdf.read_generic_hdf5(os.path.join(fpath))
      data_vol.append(data)

      #Conversao para refletividade
      data_raw = data["dataset1/data1/data"]["data"]
      gain = data["dataset1/data1/what"]["attrs"]["gain"]
      offset = data["dataset1/data1/what"]["attrs"]["offset"]
      dbZ = data_raw * gain + offset

      for i in range(len(dbZ)):
        for j in range(len(dbZ[i])):
          if dbZ[i][j] <= 10:
            dbZ[i][j] = 0

      #th_vals.append(attrs['Elevation'])
      print(data.keys())
      #print(data["dataset3/data1/data"].keys())
      #print(data["dataset3/data1/data"]["attrs"])
      #print(data["dataset3/data1/data"]["data"].shape)
      fig = plt.figure(figsize=(10,10))
      da = wrl.georef.create_xarray_dataarray(dbZ).wrl.georef.georeference()
      im = da.wrl.vis.plot(fig=fig, crs="cg")
  data_vol = np.array(data_vol)
  #th_vals = np.array(th_vals)

#Acessa a pasta que contem os dados e cria um gif com os plots do radar
def gerarGif(datalist):
  arquivos = timeFilter(typeFilter(datalist))
  for arquivo in arquivos:
    fpath = '/content/drive/MyDrive/data/' + arquivo #Alterar o caminho da pasta, caso nao seja executado no colab
    frames = []
    data = wrl.io.hdf.read_generic_hdf5(os.path.join(fpath))

    for i in range(1, 14):
      fig = plt.figure(figsize=(10,10))
      da = wrl.georef.create_xarray_dataarray(data["dataset"+str(i)+"/data1/data"]["data"]).wrl.georef.georeference()
      im = da.wrl.vis.plot(fig=fig, crs="cg")

      aux = "dataset"+str(i)
      frame_path = f'/content/drive/MyDrive/frames/{aux}.png' #Alterar o caminho da pasta, caso nao seja executado no colab
      plt.savefig(frame_path, dpi=100)
      frames.append(imageio.imread(frame_path))
      plt.close()

      gif_path = '/content/drive/MyDrive/gifs/'+arquivo+'.gif' #Alterar o caminho da pasta, caso nao seja executado no colab
      imageio.mimsave(gif_path, frames, fps=1, loop=0)

#Codigo para determinar a precipitacao a partir dos dados obtidos
def precipitacao(datalist):
  #Coordenadas
  radar = (-20.27855, -54.47396) #Latitude e Longitude do radar em Jaraguari
  estacao = (-20.45, -54.72) #Latitude e Longitude da estacao pluviometrica em Campo Grande

  #Distancia (em metros)
  distancia_m = geodesic(radar, estacao).meters
  #print(f"Distancia: {distancia_m:.2f} m")

  #Azimute
  g = Geodesic.WGS84.Inverse(radar[0], radar[1], estacao[0], estacao[1])
  azimute = g['azi1'] % 360
  #print(f"Azimute: {azimute:.2f}º")

  #Lista de taxas de precipitacao
  R_Norte = []
  R_Sudeste = []

  for arquivo in datalist:
      fpath = '/content/drive/MyDrive/data/' + arquivo #Alterar o caminho da pasta, caso nao seja executado no colab

      #Extracao dos dados
      data = wrl.io.hdf.read_generic_hdf5(os.path.join(fpath))

      for i in range(1, 14):
          #Conversao para refletividade
          data_raw = data["dataset" + str(i) + "/data1/data"]["data"]
          gain = data["dataset" + str(i) + "/data1/what"]["attrs"]["gain"]
          offset = data["dataset" + str(i) + "/data1/what"]["attrs"]["offset"]
          dbZ = data_raw * gain + offset
          #print(data_raw)
          #print(dbZ)

          #Plot do radar
          #fig = plt.figure(figsize=(10,10))
          #da = wrl.georef.create_xarray_dataarray(dbZ).wrl.georef.georeference()
          #im = da.wrl.vis.plot(fig=fig, crs="cg")
          #print(data.keys())
          #print(data["where"]["attrs"]["height"])

          #Parametros do radar
          rscale = data["dataset" + str(i) + "/where"]["attrs"]["rscale"]
          nrays = data["dataset" + str(i) + "/where"]["attrs"]["nrays"]
          nbins = data["dataset" + str(i) + "/where"]["attrs"]["nbins"]
          endtime = data["dataset" + str(i) + "/what"]["attrs"]["endtime"]

          #Indices na matriz
          az_idx = int(round(azimute)) % nrays
          bin_idx = int(round(distancia_m / rscale))

          valor_dbz = dbZ[az_idx, bin_idx]
          if valor_dbz <= 0:
              valor_dbz = 0.0
              Z = 0.0
          else:
              Z = 10**(valor_dbz / 10)

          R1 = (Z / 200) ** (1/1.6) #Relacao Z-R para o estado de Sao Paulo
          R2 = (Z / 315) ** (1/1.2) #Relacao Z-R para o estado do Pará

          R_Sudeste.append(float(R1))
          R_Norte.append(float(R2))

          #print(f"Indice de azimute: {az_idx}")
          #print(f"Indice de alcance: {bin_idx}")
          print(f"Horário Leitura: {endtime}, Valor do dbZ: {valor_dbz}, R1: {R1:.2f} mm/h, R2: {R2:.2f} mm/h")

  #Transforma o acumulado de mm/h para mm
  acumulado_Sudeste = sum([r * (5/60) for r in R_Sudeste])
  acumulado_Norte = sum([r * (5/60) for r in R_Norte])

  print("\n")
  print(f"Acumulado utilizando a relação Z-R para o estado de São Paulo: {acumulado_Sudeste}")
  print(f"Acumulado utilizando a relação Z-R para o estado do Pará: {acumulado_Norte}")

#Acessa a pasta data e cria uma lista com os arquivos que estao la
datalist = os.listdir('/content/drive/MyDrive/data') #Alterar o caminho da pasta, caso nao seja executado no colab
print(datalist)
plotRadar(datalist)
#gerarGif(datalist)
#precipitacao(datalist)