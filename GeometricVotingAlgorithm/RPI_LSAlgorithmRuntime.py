'''
LIS_Algorithm_Runtime.py
Based on: Geometric Voting Algorithm for Star Trackers
Developed by Ignacio Bugueno
Contact: ignacio.bugueno@live.com
Supervised by Samuel Gutierrez
'''

import pandas as pd
import math
import os
import subprocess
from PIL import Image
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from collections import Counter

#Load catalog
df_catalog = pd.read_csv('tmp/Catalog_DataFrame.csv', sep = '\t')
df_catalog = df_catalog.drop('Unnamed: 0', 1)

#Load distance table
df_T = pd.read_csv('tmp/T_DistanceTable_DataFrame.csv', sep = '\t')
df_T = df_T.drop('Unnamed: 0', 1)

##LIS algorithm - Running

###Source extractor

Cur_Dir = os.path.dirname(os.path.realpath('__file__')) + '/'

# Img .jpg name.
#img_jpg_name = 'img.jpg'
nombre_img_jpg = '26_07_-_21_01_51_image7_800.jpg'
# Directorio imagen .jpg.
dir_img_jpg = './Sample_images/' + nombre_img_jpg
# Directorio donde se guarda imagen .fits.
dir_img_fits = Cur_Dir
# Nombre imagen .fits.
nombre_img_fits = 'img_fits.fits'
# Directorio SExtractor.
dir_sext = './sextractor'
# Directorio base del catalogo proyectado.
path_base = './Catalog/Projected/'
# Directorio donde se guarda el archivo 'sext'.
path_stars = dir_img_fits + 'sext'
# Directorio base de catalogo no proyectado.
cat_no_proy = './Catalog/Normal/'
# Directorio donde se guarda el catalogo proyectado en el punto que entrega como resultado el primer Match.
new_path_catalog = dir_img_fits + 'new_cat'

current_path = os.path.dirname(os.path.realpath('__file__'))


## 4.- Transformacion imagen a .fits.

# Lectura de imagen.
image = Image.open(dir_img_jpg)
# Conversion a escala de grises.
imagebn = image.convert('L')
# Obtiene el tamano de la imagen.
xsize, ysize = imagebn.size
# Toma los datos de cuentas de la imagen (0 - 255).
fits_aux1 = imagebn.getdata()
# Guarda esas cuentas en un arreglo.
fits_aux2 = np.array(fits_aux1)
# Transforma ese arreglo a las mismas dimensiones de la imagen.
fits_aux3 = fits_aux2.reshape(ysize, xsize)
# Invierte el array para quedar en la orientacion adecuada.
fits_aux4 = np.flipud(fits_aux3)
# Crea un archivo .fits basico.
fits_aux5 = fits.PrimaryHDU(data=fits_aux4)
# Guarda el archivo en formato "fits".
fits_aux5.writeto(nombre_img_fits, clobber=True)

## 5.- Ejecuta SExtractor.

# Cambia el directorio al de SExtractor.
os.chdir(dir_sext)
# Define el directorio de la imagen.
imdir = dir_img_fits + nombre_img_fits
# Define el sextractor.
sext = 'sextractor ' + imdir
# Se corre sextractor. Genera "test.cat" y se lee como tabla.
subprocess.check_output(sext, shell=True)
sex_aux1 = ascii.read('./test.cat', format='sextractor')
# Ordena por magnitud y selecciona las 40 estrellas mas brillantes.
sex_aux1.sort(['MAG_ISO'])
sex_aux2 = sex_aux1[0:40]
# Cambia el directorio.
os.chdir(dir_img_fits)
## Define posiciones X e Y a partir de los datos entregados por SExtractor.
sex_x = sex_aux2['X_IMAGE']
sex_y = sex_aux2['Y_IMAGE']
sex_mag = sex_aux2['MAG_ISO']
# Conversion a coordenadas de CMOS.
sex_x1 = (sex_x - 512)*0.00270   # queda en mm centrada en (0,0)
sex_y1 = (sex_y - 512)*0.00270   # queda en mm centrada en (0,0)

os.chdir(current_path)

# Guarda las columnas X, Y y MAG del resultado de Sextractor.
ascii.write([sex_x1, sex_y1, sex_mag], 'tmp/sext', delimiter = '\t', format = 'no_header', formats = {'col0':'% 15.10f', 'col1':'% 15.10f', 'col2':'% 15.10f'})


###Read Source extractor 
with open('tmp/sext') as f:
    content = f.readlines()
    
content = [x.strip() for x in content] 

x_sext_list = []
y_sext_list = []
ID_sext_list = []

print 'First for'

for i in range(0, 15): #len(content)):
    print 'i: ' + str(i)
    sext_x, sext_y, sext_mag = content[i].split('\t')

    #Mas fuerte a magnitud 6
    if (float(sext_mag) < 6):
    
        x_sext_list.append(int(float(sext_x)/0.00270))
        y_sext_list.append(int(float(sext_y)/0.00270))
        ID_sext_list.append(i + 1)
    
df_sext = pd.DataFrame({'x': x_sext_list,
                        'y': y_sext_list,
                        'ID': ID_sext_list})

numberPossibleStars = df_sext.size/3

V_Voting_List = []

print 'Second for'
print 'numberPossibleStars: ' + str(numberPossibleStars)

for i in range (0, numberPossibleStars):
    print 'i: ' + str(i)
    V_i = []
    for j in range (i + 1, numberPossibleStars - 1):
        
        x_Star_i = df_sext.iloc[i].x
        x_Star_j = df_sext.iloc[j].x
        y_Star_i = df_sext.iloc[i].y
        y_Star_j = df_sext.iloc[j].y
        
        #ERR_degrees_Star_i = 
        #ERR_degrees_Star_j = 
        #e_ij_degrees = abs(ERR_degrees_Star_i + ERR_degrees_Star_j)
        e_ij = 0.25 #0.25 equivalentr a a 5 pixeles error
        
        d_ij = math.radians((48.8/1024)*math.sqrt((x_Star_i - x_Star_j)**2 + (y_Star_i - y_Star_j)**2))

        if (math.degrees(d_ij) <= 50): #62 Raspberry Pi Camera Angle of View: horizontal
            
            Rij_min = d_ij - e_ij
            Rij_max = d_ij + e_ij
            
            for star_i in df_T['ID1'].where((df_T['d_radians'] >= Rij_min) & (df_T['d_radians'] <= Rij_max)):
                if (math.isnan(star_i) == False):    
                    V_i.append(star_i)
            
            for star_j in df_T['ID2'].where((df_T['d_radians'] >= Rij_min) & (df_T['d_radians'] <= Rij_max)):       
                if (math.isnan(star_j) == False):    
                    V_i.append(star_j)
        
    V_Voting_List.append(V_i)    

#for all possible stars Si do

def most_common(lst):
    if (lst == []):
        return -1
    else:    
        return max(set(lst), key=lst.count)

St_i_list = []

print 'Third for'
print 'len(V_Voting_List): ' + str(len(V_Voting_List))

i = 0

for Si_sublist in V_Voting_List:
    print 'i: ' + str(i)
    if (Si_sublist == []):
        St_i_list.append(-1)
    else:    
        c = Counter(Si_sublist)
        St_i_list.append(c.most_common(1)[0][0])
        #St_i_list.append(most_common(Si_sublist)) #Assign St_i to the catalogue star which got the maximal number of votes
    i = i + 1

print St_i_list

St_i_list = map(int, St_i_list)

Match_i = []
Match_j = []

numberMatchs = 0

print 'Fourth for'
print 'numberPossibleStars: ' + str(numberPossibleStars)

for i in range (0, numberPossibleStars):
    print 'i: ' + str(i)
    for j in range (i + 1, numberPossibleStars - 1):
        
        if (St_i_list[i] != -1 and St_i_list[j] != -1):

            DEC_radians_Star_i = math.radians(df_catalog.iloc[St_i_list[i] - 1].DEC_degrees)
            DEC_radians_Star_j = math.radians(df_catalog.iloc[St_i_list[j] - 1].DEC_degrees)
            
            RA_radians_Star_i = math.radians(df_catalog.iloc[St_i_list[i] - 1].RA_degrees)
            RA_radians_Star_j = math.radians(df_catalog.iloc[St_i_list[j] - 1].RA_degrees)
            
            e_ij_degrees = 0.25

            d_i = DEC_radians_Star_i
            d_j = DEC_radians_Star_j
            a_i = RA_radians_Star_i
            a_j = RA_radians_Star_j

            d_ij = abs(math.acos(math.sin(d_i) * math.sin(d_j) + math.cos(d_i) * math.cos(d_j) * math.cos(a_i - a_j))) #Angular distance Dij = |Si - Sj|

            Rij_min = d_ij - math.radians(e_ij_degrees)
            Rij_max = d_ij + math.radians(e_ij_degrees)

            if ((d_ij >= Rij_min) & (d_ij <= Rij_max)):
                list_i = [df_sext.iloc[i].ID, int(df_catalog.iloc[St_i_list[i] - 1].ID)]
                Match_i.append(list_i)

                list_j = [df_sext.iloc[j].ID, int(df_catalog.iloc[St_i_list[j] - 1].ID)]
                Match_j.append(list_j)

                numberMatchs = numberMatchs + 1

print 'Numero de matchs: ' + str(numberMatchs)

# Python code to remove duplicate elements
def Remove(duplicate):
    final_list = []
    for num in duplicate:
        if num not in final_list:
            final_list.append(num)
    return final_list
     
Match_i = Remove(Match_i)
Match_j = Remove(Match_j)

Match_list = Match_i + Match_j

Match_list = Remove(Match_list)

for pair in Match_list:
    print 'Match: Star ' + str(pair[0]) 
    print 'Match ID in Catalog: ' + str(pair[1])
    print 'ID: ' + str(df_catalog.iloc[pair[1] - 1].ID)
    print 'RA: ' + str(df_catalog.iloc[pair[1] - 1].RA_degrees)
    print 'DEC: ' + str(df_catalog.iloc[pair[1] - 1].DEC_degrees)
    print ''

#Estimate attitude using a least squares quaternion based method

