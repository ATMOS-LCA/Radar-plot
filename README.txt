O código inicial para realizar o plot dos dados fornecidos pelo radar de Jaraguari possui 5 funções: typeFilter(), timeFilter(), gerarGif(), plotRadar() e precipitação().

Os arquivos utilizados para testar o programa foram retirados do site da CEMADEN, no formato .vol.h5. Esses arquivos são nomeados da seguinte forma: yyyymmddhhmmss00tipo.vol.h5, no qual:
- yyyymmdd é o dia que foi realizado a leitura;
- hhmmss é a hora que foi realizado a leitura;
- tipo é o tipo de dado armazenado no arquivo;
Exemplo: 2025031818400300dBZ.vol.h5

Função typeFilter(datalist):
	Essa função recebe a lista com os nomes dos arquivos retirados do site da CEMADEN e filtra os arquivos pelo tipo de dado que se deseja analisar. Pelo exemplo acima, é possível verificar que o tipo de dado armazenado no arquivo é dBZ.

Função timeFilter(datalist):
	Essa função recebe a lista com os nomes dos arquivos retirados do site da CEMADEN e filtra os arquivos pelo período de tempo especificado. Como mostrado anteriormente, podemos saber a data da leitura realizada pelo radar através do nome do arquivo.

Função plotRadar(datalist):
	Essa função recebe a lista com os nomes dos arquivos retirados do site da CEMADEN e gera um gráfico dos dados contidos nos arquivos. Nessa função, é possível realizar o plot filtrando os arquivos que se deseja analisar, utilizando as funções typeFilter() e timeFilter().

Função gerarGif(datalist):
	Essa função recebe a lista com os nomes dos arquivos retirados do site da CEMADEN e gera um gif, combinando os plot do dados armazenados nos arquivos. Também é possível realizar o plot filtrando os arquivos, utilizando as funções typeFilter() e timeFilter().

Função precipitação(datalist):
	Essa função recebe a lista com os nomes dos arquivos retirados do site da CEMADEN e determina a precipitação ocorrida a partir dos dados obtidos. Para isso, foi necessário localizar tanto o radar de Jaraguari quanto a estação pluviométrica de Campo Grande na matriz de dados de dBZ, converter o dado em dBZ para Z e aplicar a função inversa da relação Z-R que se deseja testar. Ao final, é realizado uma soma desses valores, multiplicando cada valor pelo tempo de leitura do radar.
