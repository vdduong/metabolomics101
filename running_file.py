# running_file.py
# running the algorithm

import matplotlib.pyplot as plt
from functions import *
import datetime


print '*'*10
### here are the values of tolerance for matching and alignment
tol_H = 0.08
tol_alignment = 0.05
tol_hausdorff = tol_H

name_file_spectrum = "/Users/James/Metabo_1/Exp_peak_list/43-sn10.peaks"
peak_list = spectrum(name_file_spectrum)

# for the peak list plotting
#list_peak =[]

list_file = data_file()

###

### for the plotting of theoretical and experimental patterns of metabolites

###
# for the peak list plotting
#list_peak =[]
#list_file = list()
###
#name_file_1 = "/Users/James/Metabo_1/database/bmse000408.str"
#list_file.append(name_file_1)


### Manual implementation of database to boost the speed :
#list_file.append("/Users/James/Metabo_1/database/bmse000931.str") # argininosuccinate
#list_file.append("/Users/James/Metabo_1/database/bmse000155.str") # creatinine
#list_file.append("/Users/James/Metabo_1/database/bmse000610.str") # N-acetylglycine
#list_file.append("/Users/James/Metabo_1/database/bmse000977.str") # glycine
#list_file.append("/Users/James/Metabo_1/database/bmse000855.str") # D-Glucose
#list_file.append("/Users/James/Metabo_1/database/bmse000076.str") # citrate
#list_file.append("/Users/James/Metabo_1/database/bmse000038.str") # glutamine 
#list_file.append("/Users/James/Metabo_1/database/bmse000811.str") # valine
#list_file.append("/Users/James/Metabo_1/database/bmse000920.str") # leucine
#list_file.append("/Users/James/Metabo_1/database/bmse000866.str") # isoleucine
#list_file.append("/Users/James/Metabo_1/database/bmse000361.str") # 2-hydroxybutyrate
#list_file.append("/Users/James/Metabo_1/database/bmse000426.str") # TMAO
#list_file.append("/Users/James/Metabo_1/database/bmse000948.str") # betaine
#list_file.append("/Users/James/Metabo_1/database/bmse000398.str") # methylmalonate
#list_file.append("/Users/James/Metabo_1/database/bmse000152.str") # trigonelline
#list_file.append("/Users/James/Metabo_1/database/bmse000994.str") # alanine
#list_file.append("/Users/James/Metabo_1/database/bmse000120.str") # taurine
#list_file.append("/Users/James/Metabo_1/database/bmse000915.str") # methionine
#list_file.append("/Users/James/Metabo_1/database/bmse000900.str") # phenylalanine
#list_file.append("/Users/James/Metabo_1/database/bmse000051.str") # tyrosine
#list_file.append("/Users/James/Metabo_1/database/bmse000818.str") # lactate
#list_file.append("/Users/James/Metabo_1/database/bmse000715.str") # 2-hydroxyphenylacetate
#list_file.append("/Users/James/Metabo_1/database/bmse000406.str") # glutaric acid


#m= matrix_adjacent(name_file_1)
#shift_table =  shift(name_file_1)
#pattern_tocsy = pattern(shift_table,m)
#dict_assignment = assignment(pattern_tocsy, peak_list, tol_H)
#for key in dict_assignment.keys():
#  print key,dict_assignment[key]
nb_components = 0
list_point, list_point_1 = [], []
for name_file in list_file :
    nb_components, list_point, list_point_1 = \
               iterative_matching(name_file, peak_list, tol_H, tol_alignment, \
                                  nb_components, list_point, list_point_1)
print nb_components

