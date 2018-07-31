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

#data={'taxonomy':[],'imt':[],'noDamageLimit':[],'iml':[],'d1':[],'d2':[],'d3':[],'d4':[],'d5':[]}
data={'taxonomy':[],'imt':[],'iml_unit':[],'noDamageLimit':[],'iml':[],'d0':[],'d1':[],'d2':[],'d3':[],'d4':[],'d5':[]}
#fieldnames=['taxonomy','imt','iml_unit','noDamageLimit','iml','d0','d1','d2','d3','d4','d5']

#given the ground motion for which the poes should be calculated III-XI in 0.5 steps
gm = [3.+i*.5 for i in range(17)]

#read taxonomy and vi from file
with open('frag_vis.csv','r') as csvfile:
    headerlines = 1
    delimiter = ','
    reader = csv.reader(csvfile,delimiter=delimiter)
    #skip header
    for i in range(headerlines):
        next(reader,None)

    for row in reader:
        #get one row
        #remove whitespace
        row = [x.strip() for x in row]

        bt_VI = float(row[1])

        #calculate mu_d for each gm step
        mu_d=[2.5*(1+math.tanh((gm_step + 6.25*bt_VI - 13.1)/2.3)) for gm_step in gm]

        #calculate damage grade distribution
        def pbeta(mu_d,t=8):
            #bounds of dg0-5
            a,b = 0,6
            #r,t parameters according to Giovinazzi 2005 (t=8)
            r = [t * (0.007 * mu_di**3 - 0.0525 * mu_di**2 + 0.2875 * mu_di) for mu_di in mu_d]
            #convert r,t to alpha,beta
            alpha = r
            beta = [t - al for al in alpha]
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

        #store everything
        data['taxonomy']=data['taxonomy']+ [row[0] for i in range(len(gm))]
        data['imt']=data['imt']+ ['MMI' for i in range(len(gm))]
        data['iml_unit']=data['iml_unit']+ ['MMI' for i in range(len(gm))]
        data['noDamageLimit']=data['noDamageLimit']+ [1 for i in range(len(gm))]
        data['iml']=data['iml']+gm
        data['d0']=data['d0']+fragility_curves[0]
        data['d1']=data['d1']+fragility_curves[1]
        data['d2']=data['d2']+fragility_curves[2]
        data['d3']=data['d3']+fragility_curves[3]
        data['d4']=data['d4']+fragility_curves[4]
        data['d5']=data['d5']+fragility_curves[5]


#get coburn and spence vulnerability function ( loss ratio fatalities for building type)
#coeff = 0.3 day 0.5 night
#fat = [round(coeff*(0.25*bt_occ[i]*dg4[i]+bt_occ[i]*dg5[i])) for i in range(len(bt))]

##probabilities for dgs 4 and 5 for each gm step
#dg4 = [dg.cdf(4.5)-dg.cdf(3.5) for dg in dg_betas]
#dg5 = [1 -dg.cdf(4.5) for dg in dg_betas]
#
##loss ratios for each gm step
#vulnerability_curves = {
#        'iml':gm,
#        'day':[],
#        'night':[]
#        }
#vulnerability_curves['day']=[0.3*(0.25*dg4[i] + dg5[i]) for i in range(len(dg4))]
#vulnerability_curves['night']=[0.5*(0.25*dg4[i] + dg5[i]) for i in range(len(dg4))]

def write_csv(data,filename):
    '''
    Write a dictionary to csv where keys of the dictionary are the column names
    '''
    with open(filename, 'w') as f:
        fieldnames=['taxonomy','imt','iml_unit','noDamageLimit','iml','d0','d1','d2','d3','d4','d5']
        #fieldnames=['taxonomy','imt','noDamageLimit','iml','d1','d2','d3','d4','d5']
        writer = csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(data[fieldnames[0]])):
            row = {}
            for key in data:
                try:
                    row[key] = data[key][i]
                except:
                    print key,i
            writer.writerow(row)

#write to csv's
write_csv(data,'frag_model.csv')
#write_csv(vulnerability_curves,filename+'_vc_lr.csv')

#convert to xml
import rmtk.parsers.fragility_model_converter as fmc

fmc.csv_to_xml('frag_model.csv','frag_metadata.csv', 'frag_model.xml')

