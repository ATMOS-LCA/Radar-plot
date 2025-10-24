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
#!pip install geopandas fiona

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
import geopandas as gpd

from shapely.geometry import Polygon, mapping

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
      fpath = '/content/drive/MyDrive/IC-Vinicius/data/' + fname #Alterar o caminho da pasta, caso nao seja executado no colab
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
    fpath = '/content/drive/MyDrive/IC-Vinicius/data/' + arquivo #Alterar o caminho da pasta, caso nao seja executado no colab
    frames = []
    data = wrl.io.hdf.read_generic_hdf5(os.path.join(fpath))

    for i in range(1, 14):
      fig = plt.figure(figsize=(10,10))
      da = wrl.georef.create_xarray_dataarray(data["dataset"+str(i)+"/data1/data"]["data"]).wrl.georef.georeference()
      im = da.wrl.vis.plot(fig=fig, crs="cg")

      aux = "dataset"+str(i)
      frame_path = f'/content/drive/MyDrive/IC-Vinicius/frames/{aux}.png' #Alterar o caminho da pasta, caso nao seja executado no colab
      plt.savefig(frame_path, dpi=100)
      frames.append(imageio.imread(frame_path))
      plt.close()

      gif_path = '/content/drive/MyDrive/IC-Vinicius/gifs/'+arquivo+'.gif' #Alterar o caminho da pasta, caso nao seja executado no colab
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
      fpath = '/content/drive/MyDrive/IC-Vinicius/data/' + arquivo #Alterar o caminho da pasta, caso nao seja executado no colab

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
          if valor_dbz <= 10:
              dbZ[az_idx, bin_idx] = 0.0
              Z = 0.0
          else:
              Z = 10**(valor_dbz / 10)

          R1 = (Z / 200) ** (1/1.6) #Relacao Z-R para o estado de Sao Paulo
          R2 = (Z / 315) ** (1/1.2) #Relacao Z-R para o estado do Pará

          R_Sudeste.append(float(R1))
          R_Norte.append(float(R2))

          #print(f"Indice de azimute: {az_idx}")
          #print(f"Indice de alcance: {bin_idx}")
          #print(f"Horário Leitura: {endtime}, Valor do dbZ: {valor_dbz}, R1: {R1:.2f} mm/h, R2: {R2:.2f} mm/h")

  #Transforma o acumulado de mm/h para mm
  acumulado_Sudeste = sum([r * (5/60) for r in R_Sudeste])
  acumulado_Norte = sum([r * (5/60) for r in R_Norte])

  #print("\n")
  #print(f"Acumulado utilizando a relação Z-R para o estado de São Paulo: {acumulado_Sudeste}")
  #print(f"Acumulado utilizando a relação Z-R para o estado do Pará: {acumulado_Norte}")

  return acumulado_Sudeste, acumulado_Norte

def graficoChuva(datalist):
  # --- 1. Preparar os Dados ---

  # Gerar datas e horas (a cada hora por um dia)
  datas_horas = pd.date_range(start='2025-03-13 00:00', end='2025-03-14 00:00', freq='h')

  # Preparando os dados de chuva acumulado por horário
  acumulado_pluviometro = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.6, 0.0, 0.0, 0.0, 0.0] #Dados de chuva retirados do INMET
  acumulado_radar_R1 = [0.0]
  acumulado_radar_R2 = [0.0]

  for i in range(0, 24):
    datalist2 = []
    R1 = 0.0
    R2 = 0.0
    for arquivo in datalist:
      if i < 10:
        if arquivo[8:10] == "0"+str(i):
          datalist2.append(arquivo)
      else:
        if arquivo[8:10] == str(i):
          datalist2.append(arquivo)
    print(datalist2)
    if len(datalist2) > 0:
      R1, R2 = precipitacao(datalist2)
    acumulado_radar_R1.append(R1)
    acumulado_radar_R2.append(R2)

  # Criar um DataFrame para organizar os dados
  df_precipitacao = pd.DataFrame({
      'DataHora': datas_horas,
      'Pluviometro': acumulado_pluviometro,
      'Radar_R1': acumulado_radar_R1,
      'Radar_R2': acumulado_radar_R2
  })

  # --- 2. Criar o Gráfico com Plotly ---

  fig = go.Figure()

  # Adicionar a série do Pluviômetro
  fig.add_trace(go.Scatter(
      x=df_precipitacao['DataHora'],
      y=df_precipitacao['Pluviometro'],
      mode='lines+markers', # Mostra linhas e marcadores para cada ponto
      name='Pluviômetro (Observado)',
      line=dict(color='blue', width=2),
      marker=dict(size=6, symbol='circle')
  ))

  # Adicionar a série do Radar (Relação Z-R1)
  fig.add_trace(go.Scatter(
      x=df_precipitacao['DataHora'],
      y=df_precipitacao['Radar_R1'],
      mode='lines+markers',
      name='Radar (Estimativa R1)',
      line=dict(color='red', width=2, dash='dash'), # Linha tracejada para diferenciar
      marker=dict(size=6, symbol='square')
  ))

  # Adicionar a série do Radar (Relação Z-R2)
  fig.add_trace(go.Scatter(
      x=df_precipitacao['DataHora'],
      y=df_precipitacao['Radar_R2'],
      mode='lines+markers',
      name='Radar (Estimativa R2)',
      line=dict(color='green', width=2, dash='dot'), # Linha pontilhada para diferenciar
      marker=dict(size=6, symbol='diamond')
  ))

  # --- 3. Personalizar o Layout do Gráfico ---

  fig.update_layout(
      title='Comparação de Acumulado Horário de Precipitação: Radar vs. Pluviômetro',
      xaxis_title='Data e Hora',
      yaxis_title='Acumulado de Precipitação (mm)',
      hovermode='x unified', # Permite ver todos os valores para um dado ponto no tempo
      template='plotly_white', # Um tema limpo para o gráfico
      legend_title_text='Fonte dos Dados',
      font=dict(
          family="Arial, sans-serif",
          size=12,
          color="black"
      )
  )

  # --- 4. Exibir o Gráfico ---
  fig.show()

def plotRadar_on_MS(datalist, dbz_threshold=10):
    """
    Gera mapas de refletividade do radar (dBZ) com o mapa de Mato Grosso do Sul
    e limites municipais como fundo, usando Cartopy e Geopandas.

    Args:
        datalist (list): Lista de nomes dos arquivos .vol.h5 a serem processados.
        dbz_threshold (int): Limite inferior de dBZ para plotagem (filtro de ruído).
    """
    SHAPEFILE_PATH = '/content/drive/MyDrive/IC-Vinicius/MS_Municipios_2024.shx'

    try:
        limites_ms = gpd.read_file(SHAPEFILE_PATH)

        # 1. Checa se o CRS está ausente ou é None
        if limites_ms.crs is None:
            print("AVISO: CRS do Shapefile ausente. Assumindo WGS84 (EPSG:4326).")
        limites_ms.crs = "EPSG:4326"

        # 2. Checa se o CRS é diferente de WGS84 e o converte
        if limites_ms.crs.to_string() != 'EPSG:4326':
            print(f"Convertendo CRS de {limites_ms.crs.to_string()} para EPSG:4326...")
            limites_ms = limites_ms.to_crs(epsg=4326)

        print(f"Limites de {len(limites_ms)} municípios do MS carregados com sucesso.")

    except Exception as e:
        # Este bloco captura erros na leitura do arquivo ou na manipulação do GeoDataFrame
        print(f"ERRO CRÍTICO ao carregar o Shapefile dos limites: {e}")
        print("O mapa será plotado, mas sem os limites municipais.")
        # Garante que limites_ms seja None se houver falha, para evitar o plot
        limites_ms = None

    # Coordenadas do Radar de Jaraguari-MS (para centralizar o mapa)
    RADAR_LON = -54.47396
    RADAR_LAT = -20.27855

    # Coordenadas de Extensão do Mapa
    MIN_LON = -58.5
    MAX_LON = -50.5
    MIN_LAT = -24.5
    MAX_LAT = -17.0

    for fname in datalist:
        fpath = os.path.join('/content/drive/MyDrive/IC-Vinicius/data', fname)

        try:
            # Leitura do arquivo .vol.h5
            data = wrl.io.hdf.read_generic_hdf5(fpath)
        except Exception as e:
            print(f"Erro ao ler o arquivo {fname}: {e}")
            continue

        # --- Extração e Conversão para Refletividade ---
        data_raw = data["dataset1/data1/data"]["data"]
        gain = data["dataset1/data1/what"]["attrs"]["gain"]
        offset = data["dataset1/data1/what"]["attrs"]["offset"]
        dbZ = data_raw * gain + offset

        # --- Aplicação do Filtro dBZ (Melhoria na Visualização) ---
        # Substitui os valores abaixo do limiar por NaN
        dbZ[dbZ <= dbz_threshold] = np.nan

        # --- Georreferenciamento com Wradlib e Xarray ---
        # Usa o método do wradlib para criar um DataArray do xarray com georreferenciamento
        da = wrl.georef.create_xarray_dataarray(dbZ).wrl.georef.georeference()

        # Obter os metadados do arquivo
        try:
            rscale = data["dataset1/where"]["attrs"]["rscale"] # Resolução (em metros)
            rstart = data["dataset1/where"]["attrs"]["rstart"] # Posição inicial

        except KeyError:
            print("AVISO: Metadados de rscale ou rstart não encontrados em 'where/attrs'.")
            rscale = 250.0 # 500 metros (0.5 km)
            rstart = 0.0

        # O número de bins é o tamanho da segunda dimensão do dBZ
        nbins = dbZ.shape[1]

        # Recria o vetor de range (alcance) em metros
        r = np.arange(nbins) * rscale + rstart + (rscale / 2.0)

        # Obter os arrays de coordenadas polares 1D
        az = da["azimuth"].values
        el_scalar = da["elevation"].values[0] # Assumimos a primeira elevação (0.0)

        # Obter as coordenadas do radar
        radar_alt = data["where"]["attrs"]["height"] # Altura do radar (HCA)

        # Criar os arrays 2D (azimuth e range para cada ponto)
        # O np.meshgrid combina os arrays 1D r e az em grades 2D para cada ponto de medição
        r_2d, az_2d = np.meshgrid(r, az)

        #Criar o array 2D de Elevação com o mesmo shape
        el_2d = np.full_like(r_2d, el_scalar)

        # Calcular Lat/Lon usando wradlib
        # Agora, todos os argumentos espaciais (r, az, el) são 2D
        coords = wrl.georef.spherical_to_proj(
            r_2d,
            az_2d,
            el_2d,
            (RADAR_LON, RADAR_LAT, radar_alt)
        )

        # Desempacotar o array 3D para as grades 2D
        lon_grid = coords[..., 0] # Longitude é o primeiro índice
        lat_grid = coords[..., 1] # Latitude é o segundo índice
        # alt_grid = coords[..., 2] # Altitude é o terceiro índice (não usado no plot 2D)

        try:
            # Tenta a estrutura padrão ODIM: /what/attrs/date e /what/attrs/time
            date_str = data["what"]["attrs"]["date"]
            time_str = data["what"]["attrs"]["time"]
            scan_time = f"{date_str} {time_str}"
        except KeyError:
            # Se falhar, usa o nome do arquivo
            scan_time = f"Time Missing (File: {fname})"

        # --- Configuração e Plotagem do Mapa (Cartopy) ---

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

        # Definir a Extensão do Mapa (Foco em MS/Radar)
        ax.set_extent([MIN_LON, MAX_LON, MIN_LAT, MAX_LAT],
                      crs=ccrs.PlateCarree())

        #Adicionar Feições Geográficas Básicas
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.8, edgecolor='gray', zorder=1)

        # Plotar os Limites dos Municípios
        if limites_ms is not None:
            # Garante que os limites sejam plotados
            limites_ms.plot(ax=ax,
                            edgecolor='black',
                            facecolor='none',
                            linewidth=0.5,
                            alpha=0.7,
                            transform=ccrs.PlateCarree(),
                            zorder=2)

        # Plotar os Dados de Refletividade (Raster)
        mesh = ax.pcolormesh(lon_grid, lat_grid, dbZ,
                             cmap='jet', vmin=0, vmax=60,
                             transform=ccrs.PlateCarree())

        # 5. Adicionar o Local do Radar
        ax.plot(RADAR_LON, RADAR_LAT, 'o', color='white', markersize=5,
                markeredgecolor='black', transform=ccrs.PlateCarree(), label='Radar INMET')

        # 6. Finalização do Gráfico
        cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.05)
        cbar.set_label(f'Refletividade (dBZ > {dbz_threshold} dBZ)')

        ax.set_title(f'Refletividade em MS - {fname}\nVarredura: {scan_time}')

        # Adicionar grades e rótulos
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        plt.show()

#Acessa a pasta data e cria uma lista com os arquivos que estao la
datalist = os.listdir('/content/drive/MyDrive/IC-Vinicius/data') #Alterar o caminho da pasta, caso nao seja executado no colab
#print(datalist)
#plotRadar(datalist)
#gerarGif(datalist)
#precipitacao(datalist)
#graficoChuva(datalist)
plotRadar_on_MS(datalist)