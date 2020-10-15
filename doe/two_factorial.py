# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 17:19:34 2020

@author: Gabriel
"""


import numpy as np

import scipy.stats as stats
import matplotlib.pyplot as plt



alphabet = ['A','B','C','D','E','F','G','H']

kvpairs = [(i,char) for i,char in enumerate(alphabet)]

alphadict = {char:2**i for (i,char) in kvpairs}

 
def generate_ci_matrix(k):

       
    if k > len(alphabet):
        raise Exception('k too large, Add more letters to alphabet list')
            
    
    ci_mat = np.zeros((2**k,2**k))
    ci_mat[:,0] = 1
    effect_names = ['I']
    for i in range(1,2**k):
        
        col = ci_mat[:,i]
        
        n = np.log2(i)
        if n.is_integer():
            last = i
            pattern = [-1]*i + [1]*i
            reps = len(col)/len(pattern)
            column = pattern*int(reps)
            ci_mat[:,i] = column
            effect_names.append(alphabet[int(n)])
        else:
            j = i -last
            ci_mat[:,i] = ci_mat[:,last]*ci_mat[:,j]
            effect_names.append(effect_names[j] + alphabet[int(n)] )
            
    return ci_mat,np.array(effect_names) 


def anova(y,k):    
    r,n = y.shape
    
    ci_mat,effect_names = generate_ci_matrix(k)
    C = y.sum(axis=1)@ci_mat
    effects = C/(4*n)
    SS = C**2/(8*n)
    
    
    ydd = np.sum(y,axis=(0,1))
    SST = np.sum(y**2,axis=(0,1)) - ydd**2/(r*n)
    
    SSE = SST - SS[1:].sum()
    
    
    SS_array = np.array([*SS[1:],SSE,SST])
    df_array = np.array([*np.ones_like(SS[1:]),r*(n-1),r*n-1])
    MS_array = SS_array/df_array
    MS_array[-1] = np.nan
    F0_array = MS_array/MS_array[-2]
    F0_array[-2] = np.nan
    
    
    P0 = []
    for df,F0 in zip(df_array,F0_array):
        if np.isnan(F0):
            P0.append(np.nan)
        else:    
            P0.append(stats.f.sf(F0,df,df_array[-2]))
    
    P0_array = np.array(P0)
    
    index = [*effect_names[1:],'Error','Total']

    effects_tup = (effects, effect_names)
    return index, SS_array, df_array, MS_array, F0_array, P0_array, effects_tup
    

def print_percent(res):
    index, SS_array, df_array, MS_array, F0_array, P0_array, effects_tup = res
    effects = effects_tup[0]
    percent = SS_array/SS_array[-1]
    print('\n')
    for name,effect,SS,percent in zip(index[:-2], effects[1:],
                                      SS_array[:-2], percent[:-2]):
        print(f'{name:<5}\t {effect:>10}\t {round(SS,2):>10}\t  {round(100*percent,4):>10}')


def print_anova(res):
    index, SS_array, df_array, MS_array, F0_array, P0_array, effects_tup = res
        
    print(f'\n{"Source":<6}\t{"SS":>10}\t{"df":>4}\t{"MS":>10}\t{"F0":>10}\t{"P0"}')
    for name,ss,df,ms,f0,p0 in zip(index,SS_array,df_array,
                                   MS_array,F0_array,P0_array):
        string = f'{name:<6}\t'
        string += f'{round(ss,2):>10}\t'
        string += f'{df:>4}\t'
        string += f'{round(ms,2):>10}\t'
        string += f'{round(f0,4):>10}\t'
        string += f'{p0:.5f}'
        print(string)
        
# res = anova(y,3)
# index, SS_array, df_array, MS_array, F0_array, P0_array, effects_tup = res
# effects, effect_names = effects_tup
# print_percent(res)
# print_anova(res)

# intercept = np.average(y,axis=(0,1))

# x = (-1,1,1)
# #def reg_eq(x,effects,effect_names,intercept):
# effects = effects[1:]
# effect_names = effect_names[1:]

# x_dict = {}
# for i,letter in enumerate(effect_names[-1]):
#     x_dict[letter] = x[i]

# x_array = np.ones_like(effects)
# for i,name in enumerate(effect_names):
#     for letter in name:
#         x_array[i] *= x_dict[letter]
# x_array = np.array(x_array)

# y_pred = intercept + np.sum(effects/2*x_array) 






def abc_to_tup(abc,k):
    nums=np.zeros(k)
    for char in abc:
         nums[alphabet.index(char)] = 1
    return tuple(nums)
    

def tup_to_abc(tup):
    abc = ''
    for i,num in enumerate(tup):
        if num > 0:
            abc += alphabet[i]
    k = len(tup)
    return abc,k

    
def alias(abc,generator,k):
    if abc == 'I':
        return generator
    a_num = abc_to_tup(abc, k)
    g_num = abc_to_tup(generator, k)
    
    num = np.mod(np.array(a_num)+np.array(g_num),2)
    
    abc,k = tup_to_abc(tuple(num))
    return abc

def aliases(words,generator,k):
    alias_struct = {}
    for word in words:
        alias_struct[word] = alias(word,generator,k)
    return alias_struct

def char_to_col(char,k):
    
    if char == 'I':
        column = np.ones(2**k)
        return np.array(column)
    
    i = alphadict[char]
    pattern = [-1]*i + [1]*i
    reps = 2**k/(i*2)
    column = pattern*int(reps)
    
    
    
    if len(column) < 1:
        raise Exception('Empty column, check that you have correct word length (k)')
    return np.array(column)

def abc_to_col(abc,k):
    column = np.ones(2**k)
    for char in abc:
        if char == '-':
            column *=-1
        else:
            column *= char_to_col(char, k)
    return column

def matrix_from_words(words,k,nc=0):
    columns = []
    for word in words:
        columns.append(abc_to_col(word, k))
        
    matrix = np.array(columns,dtype=float).T
    
    for i in range(nc):
        matrix = np.vstack((matrix,np.zeros(matrix.shape[1])))
        
    return matrix



def anova_unrep(y_resp,effect_string,nc):
    try:
        r,n = y_resp.shape
    except ValueError:
        r = len(y_resp)
        n = 1
        
    N = r*n
        
    nf = r - nc 
    
    k = np.log2(nf)
    
    if not k.is_integer():
        raise Exception('Unexpected dimension of y_vector, check nc')
    k = int(k)
    
    index = effect_string.split(' ') + ['Error'] + ['Total']
    ci_matrix = matrix_from_words(effect_string.split(' '),k,nc)
    
    if n ==1:
        C = y_resp @ ci_matrix
    else:
        C = np.sum(y_resp,axis=1) @ ci_matrix
    
    #effects = C/(4*n)
    SS = C**2/(8*n)
    SST = np.sum((y_resp - np.mean(y_resp))**2)
    SSE = SST - sum(SS)             
    df = np.array([1]*len(C) + [(N-1) - len(C)] + [N -1])
    df[-2] = df[-1] - sum(df[:-2])
    SS = np.array([*SS,SSE,SST])       
    
    MS = SS/df
    MS[-1] = np.nan
    F0 = MS/MS[-2]
    F0[-2] = np.nan
    
    P0 = np.empty_like(F0)
    for i,f in enumerate(F0):
        if not np.isnan(f):
            P0[i] = stats.f.sf(f,df[i],df[-2])
        else:
            P0[i] = np.nan
    return SS,df,MS,F0,P0,index



def print_anova(anova_res,width=(10,4,10,10,10),decimals=(2,0,2,4,9)):
    SS,df,MS,F0,P0,index = anova_res
    
    arrays = [SS,df,MS,F0,P0]
    
    labels = ['SS', 'df','MS','F0','P0']
    
    header = f'\n{"Source":<8}\t'
    for j,label in enumerate(labels):
        header += f'{label:<{width[j]}}\t'
    
    print(header)
    for i,name in enumerate(index):
        s = f'{name:<8}\t'
        for j,array in enumerate(arrays):
            s += f'{round(array[i],decimals[j]):<{width[j]}}\t'
        print(s)

y = np.array([[550,604],
              [669,650],
              [633,601],
              [642,635],
              [1037,1052],
              [749,868],
              [1075,1063],
              [729,860]],dtype=float)



effect_string = 'A B C AB AC BC ABC'
y_resp = np.array([32,37,30,38,6,28,9,26,24,23,27,25],dtype=float)
nc=4
k =3



# res = anova_unrep(y_resp, effect_string, nc)
# print_anova(res)

# res = anova_unrep(y_resp, 'A C AC', nc)
# print_anova(res)


# res = anova_unrep(y,effect_string, 0)
# print_anova(res)

y_resp = np.array([-3.81,2.05,13.89,1.91,-2.95,0.47,12.67,1.02])

SS,df,MS,F0,P0,index = anova_unrep(y_resp,'A B -AB', 0)


index = ['A', 'D', '-AD', 'Error', 'Total']

res = SS,df,MS,F0,P0,index
print_anova(res)


beta = [np.average(y_resp), *np.sqrt(2*SS[:3])/4]


y_pred = lambda x: beta[0] -beta[1]*x[0] + beta[2]*x[1] + beta[3]*x[2]


x = matrix_from_words('A -AB B'.split(' '), 3)


residual = []
predicted = []
for i in range(len(y_resp)):
    print(x[i])
    e = y_resp[i] - y_pred(x[i])
    residual.append(e)
    predicted.append(y_pred(x[i]))
    print(e)

plt.plot(predicted,residual,'.')




for e,row in zip(residual,x):
    print(row,e)
    
    
plt.plot(x[:,0],residual,'o')
    
    