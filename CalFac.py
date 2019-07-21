import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
tmp1=sys.stdout ###### only if one wants to shut down the useless output and warning
tmp2=sys.stderr ######
sys.stdout = open("/dev/null", "a")
sys.stderr = open("/dev/null", "a")
from astroquery.gaia import Gaia
#import warnings
#warnings.filterwarnings("ignore")

def gaia(ra,dec,radius1):
    coord = SkyCoord(ra=ra, dec=dec, unit = "deg")
    radius = u.Quantity(radius1, u.deg)
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()
    return r['phot_rp_mean_mag','phot_g_mean_mag','ra','dec']

class profile:
    def __init__(self,ra,dec):
        self.ra=ra
        self.dec=dec
                
    def Fit(self):
        data=np.loadtxt("flux_influenceVSradius.txt",dtype=float)
        x=data[:,0]
        y=data[:,1]
        yerr=data[0:,2]
        index=yerr<0.15
        x=x[index]
        y=y[index]
        yerr=yerr[index]
        total=np.column_stack((x,y,yerr))
        total=total[np.lexsort(total[:,::-1].T)]
        x=total[:,0]
        y=total[:,1]
        yerr=total[:,2]
        self.x=x
        self.y=y
        self.yerr=yerr
        yfit=np.polyfit(x,y,10)
        y3=np.poly1d(yfit)
        self.yfit=y3
    def PlotRel(self):
        fig,ax=plt.subplots(figsize=(8,5))
        plt.xlabel("Distance in Pixels",fontsize=16)
        plt.ylabel("Relative Flux",fontsize=16)
        plt.title('Relative Flux in Certain Distance',fontsize=20)
        ax.errorbar(self.x,self.y,yerr=self.yerr,color="b",fmt='o',markersize='2',elinewidth=0.5)
        plt.plot(self.x,self.yfit(self.x),"r--")
        plt.show()
    def DeBlending(self):
        job=gaia(self.ra,self.dec,0.0058333*6)
        data1=np.array([job['phot_rp_mean_mag'],job['phot_g_mean_mag'],job['ra'],job['dec']]).T
        rband=data1[:,0]
        gband=data1[:,1]
        rband=np.nan_to_num(rband)
        dif=rband-gband
        diff=np.median(dif)
        i=1
        while i<len(dif):
            if rband[i]==0:
                rband[i]=gband[i]+diff
            i+=1
        f=[]
        ferr=[]
        i=1
        #error=np.median(self.yerr)
        while i<len(dif):
            dis=np.cos(data1[i,3]*np.pi/180)*np.cos(data1[i,3]*np.pi/180)*np.power(data1[i,2]-data1[0,2],2)+np.power(data1[i,3]-data1[0,3],2)
            dis=np.sqrt(dis)*3600/21
            #index=(self.x<dis+0.5) & (self.x>dis-0.5)
            error=np.median(self.yerr)
            factor=self.yfit(dis)
            flux=np.power(10,(rband[0]-rband[i])/2.5)
            f=np.append(f,flux*factor)
            ferr=np.append(ferr,flux*error)
            i+=1
        sys.stdout=tmp1
        sys.stderr=tmp2
        totblend=np.sum(f)
        self.fact=totblend
        print ("flux fraction")
        print (totblend)
        toterror=np.sqrt(np.sum(ferr*ferr))
        print ("fraction uncertainty")
        print (toterror)
        

if __name__=='__main__':        

    #ra=111.509499
    #dec=7.615806
    #ra=32.781708
    #dec=2.418000
    ra=107.600275
    dec=-39.097393
    ################# fit the relation between flux fraction and distance###############
    
    Rel=profile(ra,dec)
    Rel.Fit()
    Rel.PlotRel()
    ######################use the fit result to calculate the blending flux fraction
    Rel.DeBlending()
    #print (Rel.fact) attribute Rel.fact is the flux fraction needed. 
    
    #######the true flux of source in TESS is following: f_true=f_observed/(1+factor), factor is the flux fraction derived in this code. 



