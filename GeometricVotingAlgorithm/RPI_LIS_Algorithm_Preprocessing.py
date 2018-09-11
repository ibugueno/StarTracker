'''
LIS_Algorithm_Preprocessing.py
Based on: Geometric Voting Algorithm for Star Trackers
Developed by Ignacio Bugueno
Contact: ignacio.bugueno@live.com
Supervised by Samuel Gutierrez
'''

import pandas as pd
import math

#Preprocessing Catalog

##Read Catalog and write Panda Dataframe

#Catalog parameters: RA - DEC - Mag

with open('catalogo/catalogo_bsc_mag_4') as f:
    content = f.readlines()
    
content = [x.strip() for x in content] 

RA_degrees_catalog_list = []
DEC_degrees_catalog_list = []
MAG_catalog_list = []
ID_catalog_list = []

for i in range(0, len(content)):
    
    RA_degrees_catalog, DEC_degrees_catalog, MAG_catalog, ID_catalog = content[i].split(" ")

    RA_degrees_catalog_list.append(float(RA_degrees_catalog))
    DEC_degrees_catalog_list.append(float(DEC_degrees_catalog))
    MAG_catalog_list.append(float(MAG_catalog))
    ID_catalog_list.append(int(ID_catalog))
    
df_catalog = pd.DataFrame({'ID': ID_catalog_list,
                           'RA_degrees': RA_degrees_catalog_list,
                           'DEC_degrees': DEC_degrees_catalog_list,
                           'MAG': MAG_catalog_list})

##Prepocessing Catalog

numberStarsInCatalog = df_catalog.size/4

T_ID1_list = []
T_ID2_list = []
T_d_list = []

for i in range (0, numberStarsInCatalog):
    print 'i: ' + str(i)
    for j in range (i + 1, numberStarsInCatalog - 1):

        DEC_radians_Star_i = math.radians(df_catalog.iloc[i].DEC_degrees)
        DEC_radians_Star_j = math.radians(df_catalog.iloc[j].DEC_degrees)
        RA_radians_Star_i = math.radians(df_catalog.iloc[i].RA_degrees)
        RA_radians_Star_j = math.radians(df_catalog.iloc[j].RA_degrees)

        d_i = DEC_radians_Star_i
        d_j = DEC_radians_Star_j
        a_i = RA_radians_Star_i
        a_j = RA_radians_Star_j

        D_ij = abs(math.acos(math.sin(d_i) * math.sin(d_j) + math.cos(d_i) * math.cos(d_j) * math.cos(a_i - a_j))) #Angular distance Dij = |Si - Sj|

        if (math.degrees(D_ij) <= 50): #62 Raspberry Pi Camera Angle of View: horizontal 
            #Append entry (i,j,Dij) to table T(ID1,ID2,d)

            T_ID1_list.append(df_catalog.iloc[i].ID)
            T_ID2_list.append(df_catalog.iloc[j].ID)
            T_d_list.append(D_ij)
            
df_T = pd.DataFrame({'ID1': T_ID1_list,
                     'ID2': T_ID2_list,
                     'd_radians': T_d_list})

df_T.sort_values('d_radians', inplace=True)

df_catalog.to_csv('tmp/Catalog_DataFrame.csv', sep='\t')
df_T.to_csv('tmp/T_DistanceTable_DataFrame.csv', sep='\t')
