# -*- coding: utf-8 -*-
"""
Created on Sun May  3 21:46:24 2020

@author: David
"""

import os
os.chdir('G:\\CPTAC')

dic_al = {}
 with open('alter_data\\maf\\TCGA.OV.mutect.b22b85eb-2ca8-4c9f-a1cd-b77caab999bd.DR-10.0.somatic.maf') as maf:
    i = 1
    while i < 7:
        i+=1    
        print(maf.readline())
    
    while True:
        x = maf.readline()
        if x:
            if x.split('\t')[0] in dic_al.keys():
                dic_al[x.split('\t')[0]].append(x.split('\t')[15][:16])
            else:
                dic_al[x.split('\t')[0]] = []
                dic_al[x.split('\t')[0]].append(x.split('\t')[15][:16])
        else:
            break
    maf.close()
id_cohort = []
for key in dic_al:
    for i in dic_al[key]:
        if i in id_cohort:
            pass
        else:
            id_cohort.append(i)
            
result = {}
for key in dic_al:
    dict = {}
    for i in dic_al[key]:
        dict[i] = dict.get(i,0)+1
    result[key] = dict    

with open('alter_MAF.csv','w') as ouf:
    ouf.write('GENESYMBOL,')
    ouf.write((',').join(id_cohort)+'\n')
    for key in result:
        ouf.write(key+',')
        for i in id_cohort[:-1]:
            if i in result[key]:
                ouf.write(str(result[key][i])+',')
            else:
                ouf.write('-'+',')
        if id_cohort[-1] in result[key]:
            ouf.write(str(result[key][id_cohort[-1]])+'\n')
        else:
            ouf.write('-'+'\n')
    ouf.close()            
    
    
    
    
        