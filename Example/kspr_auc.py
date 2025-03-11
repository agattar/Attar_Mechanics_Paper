import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import trapezoid as trapz
import os 

def data_org(file_path):

    file = open(file_path, "r")

    file = file.readlines()

    timepoints = []
    types = []
    ids = []    
    x = []
    y = []
    z = []  

    for ind, row in enumerate(file):

        if row == "ITEM: TIMESTEP\n":

            time = file[ind+1].strip()

        if len(row.split(" ")) == 5:

            timepoints.append(time)
            ids.append(row.strip().split(" ")[0])
            types.append(row.strip().split(" ")[1])
            x.append(row.strip().split(" ")[2])
            y.append(row.strip().split(" ")[3]) 
            z.append(row.strip().split(" ")[4])

    timepoints = np.array(timepoints).astype(int)
    ids = np.array(ids).astype(int)
    x = np.array(x).astype(float)
    y = np.array(y).astype(float)
    z = np.array(z).astype(float)
    types = np.array(types).astype(int) 

    df = pd.DataFrame({"time": timepoints, "id": ids, "type": types, "x": x, "y": y, "z": z})

    return df

import warnings
warnings.filterwarnings('ignore')

ksprings = np.zeros((10, 5))

auc_list = np.zeros((10, 5))

for rep in range(1, 11):

    path = f'r{rep}/dump.pulling'

    df = data_org(path)

    lamins_df = df[(df["type"] == 1) | (df["type"] == 9) | (df["type"] == 10)]

    radial_distance = (lamins_df["x"]**2 + lamins_df["y"]**2 + lamins_df['z']**2)**0.5

    lamins_df["lamin radial_distance"] = radial_distance

    lamins_df = lamins_df.sort_values(by = "id")

    lamins_df.reset_index(drop = True, inplace = True)

    chromatin_df = df[(df["type"] != 1) & (df["type"] != 9) & (df["type"] != 10)]

    radial_distance = (chromatin_df["x"]**2 + chromatin_df["y"]**2 + chromatin_df['z']**2)**0.5

    chromatin_df["chromatin radial_distance"] = radial_distance

    chromatin_df = chromatin_df.sort_values(by = "id")

    chromatin_df.reset_index(drop = True, inplace = True)

    max_y_initial = lamins_df[lamins_df['time'] == 0]['y'].max()

    max_time = lamins_df['time'].max()

    unique_times = sorted(lamins_df['time'].unique())

    strain_df = [] 
    time_df = []   
    max_y_l = []

    delta_l_l = []

    for time in unique_times:

        temp_df = lamins_df[lamins_df['time'] == time]

        max_y = temp_df['y'].max()

        strain = (max_y - max_y_initial) / max_y_initial

        delta_L = (max_y - max_y_initial)

        time_df.append(time)

        strain_df.append(strain)

        max_y_l.append(max_y)

        delta_l_l.append(delta_L)

    force_txt = f'r{rep}/Force_y.txt'

    df_force = pd.read_csv(force_txt, sep = " ", header = None)

    time_strain = pd.DataFrame({"time": time_df, "strain": strain_df, "max_y": delta_l_l})

    total = pd.concat([time_strain, df_force], axis=1)

    total.columns = ['time', 'strain', 'max_y', 'force']

    print("Replica:", rep)  

    for i in range(1,6):

        closest_value = min(time_strain['strain'], key=lambda x: abs(x - 0.1 * i))

        closest_value_PREV = min(time_strain['strain'], key=lambda x: abs(x-(i-1) * 0.1))

        print("Closest Value PREV:", closest_value_PREV)

        closest_time_prev = time_strain[time_strain['strain'] == closest_value_PREV]['time'].values[0] 

        closest_time = time_strain[time_strain['strain'] == closest_value]['time'].values[0] 

        delta_L_prev = time_strain[time_strain['strain'] == closest_value_PREV]['max_y'].values[0]

        delta_L = time_strain[time_strain['strain'] == closest_value]['max_y'].values[0]
        
        k_spring = (df_force.iloc[int(closest_time/1000)] - df_force.iloc[int(closest_time_prev/1000)]) / (delta_L - delta_L_prev)
        
        ksprings[rep-1,i-1] = k_spring
        
        index = total[total['strain'] == closest_value].index

        area_under_curve = trapz(total.iloc[:index[0],:]['force'].values,total.iloc[:index[0],:]['strain'].values)

        auc_list[rep-1,i-1] = area_under_curve

        print("Strain:", 0.1 * i)
        print("Kspring:", k_spring.values[0])
        print("Delta L:", delta_L)
        print("Closest Value:", closest_value)
        print("Closest Time:", closest_time)
        print("Area under the curve:", area_under_curve)
        print("force",df_force.iloc[int(closest_time/1000)].values[0])

        print("\n")


    file_path = "_".join(os.getcwd().split("/")[-1:] + os.getcwd().split("/")[-2:-1])

np.savetxt(f'auc_30_{file_path}.csv', auc_list, delimiter='\t')
np.savetxt(f'ksprings_30_{file_path}.csv', ksprings, delimiter='\t')

        