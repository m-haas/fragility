# Routine to create a discrete fragility function
# i.e. f(MMI,DG)=P(dg > DG|Building type)
# as well as a discrete vulnerability function
# using the Giovinazzi vulnerability index method
import math
import mcerp
import scipy.stats
import csv
import os
#and considering the vulnerability class distribution
#TODO: include uncertainty

#filename to write curves to
filename='bt5'

#given a EMS98-vulnerability class distribution for the building type
bt_vul=[0, 16 , 68 , 16 , 0, 0]
#given the ground motion for which the poes should be calculated III-X in 0.5 steps
gm = [3.+i*.5 for i in range(15)]

def vulnerabilityIndex(bt_vul):
    '''
    Returns VI according to Giovinazzi 2005 using vulnerability class shares of
    a building type
    '''
    #means&sigma of Vul Class A-F
    means = [0.90,0.74,0.58,0.42,0.26,0.10]
    sigma = 0.04
    bt_VI=0
    for vulClass in range(6):
        #normal distribution for vulClass
        VI_dist = mcerp.Normal(means[vulClass],sigma)
        #VI_dist * share, added to previous classes
        bt_VI += bt_vul[vulClass]/100*VI_dist

    return(bt_VI)

#bt_VI
bt_VI=vulnerabilityIndex(bt_vul)

#calculate mu_d for each gm step
mu_d=[2.5*(1+math.tanh((gm_step + 6.25 * bt_VI.mean - 13.1)/2.3)) for gm_step in gm]

#calculate damage grade distribution
def pbeta(mu_d,t=8):
    #bounds of dg0-5
    a,b = 0,6
    #r,t parameters according to Giovinazzi 2005 (t=8)
    r = [t * (0.007 * mu_di**3 - 0.0525 * mu_di**2 + 0.2875 * mu_di) for mu_di in mu_d]
    #convert r,t to alpha,beta
    alpha = r
    beta = [t - a for a in alpha]
    result=[]
    for i in range(len(alpha)):
        tmp=scipy.stats.beta(alpha[i],beta[i],loc=a,scale=b-a)
        result.append(tmp)

    return result

#calculate beta distributions for each gm step
dg_betas = pbeta(mu_d)

#get poes for each dg
fragility_curves = {
        'iml':gm,
        0:[],
        1:[],
        2:[],
        3:[],
        4:[],
        5:[],
        }
#pdf for each MMI level for the bt
for pdf in dg_betas:
    #for each damage grade one curve
    for dg in range(6):
        fragility_curves[dg].append(1-pdf.cdf(dg+0.5))


#get coburn and spence vulnerability function ( loss ratio fatalities for building type)
#coeff = 0.3 day 0.5 night
#fat = [round(coeff*(0.25*bt_occ[i]*dg4[i]+bt_occ[i]*dg5[i])) for i in range(len(bt))]

#probabilities for dgs 4 and 5 for each gm step
dg4 = [dg.cdf(4.5)-dg.cdf(3.5) for dg in dg_betas]
dg5 = [1 -dg.cdf(4.5) for dg in dg_betas]

#loss ratios for each gm step
vulnerability_curves = {
        'iml':gm,
        'day':[],
        'night':[]
        }
vulnerability_curves['day']=[0.3*(0.25*dg4[i] + dg5[i]) for i in range(len(dg4))]
vulnerability_curves['night']=[0.5*(0.25*dg4[i] + dg5[i]) for i in range(len(dg4))]

def write_csv(data,filename):
    '''
    Write a dictionary to csv where keys of the dictionary are the column names
    '''
    if(os.path.exists(filename)):
        print('file with that name exists none created')
    else:
        with open(filename, 'w') as f:
            fieldnames = list(data.keys())
            writer = csv.DictWriter(f,fieldnames=fieldnames)
            writer.writeheader()
            for i in range(len(data[fieldnames[0]])):
                row = {}
                for key in data:
                    row[key] = data[key][i]
                writer.writerow(row)

#write both to csv's
write_csv(fragility_curves,filename+'_fc_poes.csv')
write_csv(vulnerability_curves,filename+'_vc_lr.csv')
