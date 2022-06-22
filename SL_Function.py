import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scanpy as sc
import json
import anndata

class SL_Analysis(object):

    def __init__(self,adata,SL_pd):
        '''
        Init the analysis of synthetic lethal pair

        Parameters
        ----------
        adata:Anndata
            The raw data of scRNA-seq that calculated the pseudotime
            Tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html'
        SL_pd:pandas.DaraFrame
            DataFrame of data points with each entry in the form:['gene1','gene2']
        
        Returns
        -------
        None
        '''

        self.adata=adata
        self.SL_pd=SL_pd
        if('dpt_pseudotime' not in adata.obs.columns):
            print('......You need to calculate pseudotime at first\n')
            print('......The tutorial of pseudotime could be found at\n \
                https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html')

    def calculate_SL_count(self):
        '''
        Calculate the matrix of synthetic lethal gene

        Parameters
        ----------
        None

        Returns
        -------
        SL_count_sort:pandas.DataFrame
            the matrix of SL gene sorted by pseudotime
        '''
        print('......Analysis SL_count')
        adata=self.adata
        SL_pd=self.SL_pd
        xl=adata.var.index.tolist()
        SL=pd.DataFrame(columns=['gene1','gene2'])

        #Get the intersection of SL input and the gene name of scRNA-seq
        ret3= list(set(xl).intersection(SL_pd['gene1'].values))
        for i in ret3:
            ret4=list(set(ret3).intersection(SL_pd.loc[SL_pd['gene1']==i]['gene2'].values))
            for j in ret4: 
                SL=SL.append({'gene1':i,'gene2':j},ignore_index=True)
        ret5= list(set(xl).intersection(SL_pd['gene2'].values))
        for i in ret5:
            ret6=list(set(ret5).intersection(SL_pd.loc[SL_pd['gene2']==i]['gene1'].values))
            for j in ret6: 
                SL=SL.append({'gene1':j,'gene2':i},ignore_index=True)
        SL=SL.drop_duplicates(keep='first')

        #Drop out the same gene of SL
        sl_len=len(SL)
        for i in range(sl_len):
            if(i==sl_len):
                break
            if(SL.iloc[i]['gene1']==SL.iloc[i]['gene2']):
                print('......The same name of gene loc in',i)
                SL=SL.drop(SL.iloc[i].name)
                sl_len-=1

        #Get the result of the pair of SL
        SL_li=SL['gene1'].values.tolist()+SL['gene2'].values.tolist()
        SL_li=list(set(SL_li))
        
        #Count
        SL_count=pd.DataFrame()
        SL_count['X']=adata.obs.dpt_pseudotime
        SL_count_value=pd.DataFrame(adata[:,SL_li].X,columns=SL_li,index=adata.obs.index)
        SL_count=pd.concat([SL_count,SL_count_value],axis=1)
        self.SL_count_sort=SL_count.sort_values(by='X')
        self.SL=SL
        return self.SL_count_sort

    def calculate_SL_slope(self):
        '''
        Calculate the slope of synthetic lethal gene

        Parameters
        ----------
        None

        Returns
        -------
        slope_pd:pandas.DataFrame
            the slope result of SL gene
        '''
        #estimate SL_count_sort exist status
        try:
            self.SL_count_sort
        except NameError:
            var_exists=False
        else:
            var_exists=True
        if (var_exists==False):
            print('......You need to run calculate_SL_count first')
            return 0

        SL_count_sort=self.SL_count_sort
        print('......Analysis SL_slope')
        import scipy.stats as stats
        slope_pd=pd.DataFrame(columns=['slope', 'intercept', 'r_value', 'p_value', 'std_err','average','active','sig'])
        for i in SL_count_sort.columns[1:]: 
            if len(SL_count_sort[SL_count_sort[i]!=0].values)==0:
                slope, intercept, r_value, p_value, std_err=0,0,0,0,0
            else:  
                slope, intercept, r_value, p_value, std_err=stats.linregress(SL_count_sort[SL_count_sort[i]!=0]['X'].values, SL_count_sort[SL_count_sort[i]!=0][i].values)
            average=np.average(SL_count_sort[i].values)
            active='T'
            if(average <0.01 and average>-0.01):
                active='F'
            sig='normal'
            
            if(slope>0):
                sig='high'
            elif(slope<0):
                sig='low'
            slope_pd.loc[i]=[slope, intercept, r_value, p_value, std_err,average,active,sig]
        self.slope_pd=slope_pd
        return slope_pd
        
    
    def calculate_SL_pair(self):
        '''
        Calculate the pair of synthetic lethal gene and its type

        Parameters
        ----------
        None

        Returns
        -------
        SL_p:pandas.DataFrame
            the type of SL gene pair
        '''
        #estimate slope_pd exist status
        try:
            self.slope_pd
        except NameError:
            var_exists=False
        else:
            var_exists=True
        if (var_exists==False):
            print('......You need to run calculate_SL_slope first')
            return 0
        
        print('......Analysis SL_p')
        slope_pd=self.slope_pd
        SL=self.SL
        SL_p=pd.DataFrame(columns=['dp','type'])
        for i in range(len(SL)):
            test=[]
            g1,g2=(SL.iloc[i]['gene1']),(SL.iloc[i]['gene2'])
            #test.append(SL_count_sort[g1].values)
            #test.append(SL_count_sort[g2].values)
            #test=np.array(test)
            ty=0
            if(slope_pd.loc[g1]['active']=='T' and slope_pd.loc[g1]['sig']=='high'):
                if(slope_pd.loc[g2]['active']=='F'):
                    ty=1
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='high'):
                    ty=4
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='low'):
                    ty=3
            if(slope_pd.loc[g1]['active']=='T' and slope_pd.loc[g1]['sig']=='low'):
                if(slope_pd.loc[g2]['active']=='F'):
                    ty=2
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='high'):
                    ty=3
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='low'):
                    ty=5
            if(slope_pd.loc[g1]['active']=='F'):
                if(slope_pd.loc[g2]['active']=='F'):
                    ty=6
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='high'):
                    ty=1
                elif(slope_pd.loc[g2]['active']=='T' and slope_pd.loc[g2]['sig']=='low'):
                    ty=2

            SL_p.loc[g1+'-'+g2]={'dp':slope_pd.loc[g1]['slope']-slope_pd.loc[g2]['slope'],'type':ty}
        self.SL_p=SL_p
        return SL_p

    def save_slope_result(self,paths):
        '''
        Save the calculate result of SL

        Parameters
        ----------
        paths:list
            paths[0]:the save path of SL_count
            paths[1]:the save path of SL_slope
            paths[2]:the save path of SL_pair

        Returns
        -------
        None
        '''
        #estimate SL_p exist status
        try:
            self.SL_p
        except NameError:
            var_exists=False
        else:
            var_exists=True
        if (var_exists==False):
            print('......You need to run calculate_SL_p first')
            return 0
        self.SL_count_sort.to_csv(paths[0])
        self.slope_pd.to_csv(paths[1])
        self.SL_p.to_csv(paths[2])

    def Lazy_analysis(self,paths):
        '''
        Lazy calculate the result of SL

        Parameters
        ----------
        paths:list
            paths[0]:the save path of SL_count
            paths[1]:the save path of SL_slope
            paths[2]:the save path of SL_pair

        Returns
        -------
        None
        '''
        self.calculate_SL_count()
        self.calculate_SL_slope()
        self.calculate_SL_pair()
        self.save_slope_result(paths)

