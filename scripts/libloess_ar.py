from statsmodels.nonparametric.smoothers_lowess import lowess
from Bio.Statistics.lowess import lowess as biolowess
from scipy.interpolate import interp1d,UnivariateSpline
import numpy as np
import joblib
import scipy.integrate as integrate
import hashlib
import os
import matplotlib.pyplot as plt
import scipy
from scipy import special
from datetime import datetime

class NotImplementError(Exception):
      pass

class Myloess:
   def __init__(self,equation='cv2vslog2mu',
               use_bio_lowess=False,
               bio_lowess_f = 2.0/3.0,
               bio_lowess_iter = 3,
               use_interpolation=True,
               spline_k=4,spline_s=0.05,
               prefix=None):
       self.equation = 'cv2vslog2mu'
       self.func = None
       self.k = spline_k
       self.s = spline_s
       self.use_bio_lowess = use_bio_lowess
       self.bio_lowess_f = bio_lowess_f
       self.bio_lowess_iter = bio_lowess_iter
       self.use_interpolation = use_interpolation
       self.prefix = prefix
  
   def _convertxy(self,X,Y):
       if self.equation == "cv2vslog2mu":
           x = np.log2(np.asarray(X))
           y = np.power(np.asarray(Y)/np.asarray(X),2)
           return x,y
       elif self.equation == "cvvslog2mu":
           x = np.log2(np.asarray(X))
           y = np.asarray(Y)/np.asarray(X)
           return x,y
       else:
           raise NotImplementError("{} not implemented yet".format(self.equation))

   def findy(self,X):
       if self.equation == "cv2vslog2mu":
           newx = np.log2(np.asarray(X))
           t =np.where(np.any([newx>self.hb,newx<self.lb],axis=0))[0]
           #print("number of points out of the limit {}".format(len(t)))
           try:
               cv2 = self.func(newx)
           except ValueError:
               print("the minimum of original trainset is {:.8e}".format(np.min(self.converted_X)))
               print("the minimum of original trainset is {:.8e}".format(np.min(self.converted_Y)))
               print("current minimum = {:.8e} and current max {:.8e}".format(np.min(newx),np.max(newx)))
               raise
           #print("number of nan in cv2 is {}".format(np.sum(np.isnan(cv2))))
           cv2 = np.where(cv2>1.0e-32,cv2,1.0e-32)
           if np.min(cv2) < 0.0:
               t = np.where(cv2<0.0)[0]
               raise RuntimeError("cv2 is smaller than 0.0, \n newx = {}".format(newx[t]))
           y = np.sqrt(cv2)*X
       elif self.equation == "cvvslog2mu":
           newx = np.log2(np.asarray(X))
           t = np.where(np.any([newx>self.hb,newx<self.lb],axis=0))[0]
           cv = self.func(newx)
           cv = np.where(cv>1.0e-32,cv,1.0e-32)
           if np.min(cv) < 0.0:
               t = np.where(cv<0.0)[0]
               raise RuntimeError("cv is smaller than 0.0, \n newx = {}".format(newx[t]))
           y = cv*X
       else:
           raise NotImplementError("{} not implemented yet".format(self.equation))
       return y

   def fit_loess(self,X,Y):
       ss=self.equation + "{}".format(self.use_bio_lowess)
       if self.use_bio_lowess:
           ss += "{:02d}{:.5e}".format(self.bio_lowess_iter,self.bio_lowess_f)
       ss += "".join([" {:.10e}".format(f) for f in X])
       ss += "".join([" {:.10e}".format(f) for f in Y])
       hashstr = hashlib.md5(ss.encode()).hexdigest()
       cx,cy = self._convertxy(X,Y)
       self.X = X
       self.Y = Y
       self.converted_X = cx
       self.converted_Y = cy

       if self.prefix is not None:
           prefix = self.prefix
       else:
           prefix = "./"
       f = os.path.join("{}/{}.joblib".format(prefix,hashstr))
       if os.path.isfile(f):
           print("Myloess: {} found and used".format(f))
           res = joblib.load(f) 
           tx = res['tx']
           ty = res['ty']
       else:
           if self.use_bio_lowess:
               ty=biolowess(cx,cy,f=self.bio_lowess_f,iter=self.bio_lowess_iter)
               tx=cx
           else:
               res = lowess(cy, cx)
               tx,ty = res[:,0],res[:,1]
           res = {'tx':tx, 'ty':ty }    
           joblib.dump(res,f)
       self.loess_x = tx
       self.loess_y = ty
       if self.use_interpolation:
           self.func = interp1d(tx, ty, fill_value="extrapolate")
       else:
           self.func = UnivariateSpline(tx, ty, k=self.k, s=self.s)
       self.lb = np.min(tx)
       self.hb = np.max(tx)
      


def test_new_site(convt_object, XM, Xm):
    muM = np.mean(XM)
    mum = np.mean(Xm)
    nrep = len(muM)
    sigma = convt_object.findy(muM)/np.sqrt(nrep)
    prefactor = 1.0 / np.sqrt(2.0 * np.pi)
    new_upper_limit = (mum - muM) / sigma
    P = 2.0 * prefactor * integrate.quad(lambda x: np.exp(-np.power(x, 2) / 2.0), -np.inf, new_upper_limit)[0]
    return P

def test_new_sites(convt_object,XMs,Xms):
    muM = np.mean(XMs,axis=1)
    mum = np.mean(Xms,axis=1)
    nrep = XMs.shape[1]
    sigmas = convt_object.findy(muM)/np.sqrt(nrep)
    return sigmas,(special.erf((mum-muM)/sigmas/np.sqrt(2))+1.0)*0.5*2.0

def plot_loess(cnvt_object, prefix,labels=None):
    tx = cnvt_object.converted_X
    ty = cnvt_object.converted_Y
    print("minimum x value {:.8f}".format(np.min(tx)))
    print("maximum x value {:.8f}".format(np.max(tx)))
    xmesh = np.linspace(np.min(tx), np.max(tx), 60)
    plt.clf()     
    #plt.hexbin(tx, ty, bins=30, cmap='OrRd', gridsize=20)
    plt.hexbin(tx, ty, bins='log', cmap='OrRd', gridsize=20)
    plt.scatter(cnvt_object.loess_x,cnvt_object.loess_y,label="LOESS data")
    cb = plt.colorbar()
    cb.set_label(r'$\log_{10} N$')
    
    if labels is not None:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    plt.plot(xmesh,cnvt_object.func(xmesh), color='k',label="Fitted")
    plt.legend(loc='upper right')
    plt.savefig(prefix + "_Ndata{}.png".format(len(cnvt_object.X)))
    now = datetime.now() # current date and time
    plt.savefig("Lowess_NdataInLowess{}_TS{}.png".format(len(cnvt_object.X),now.strftime("%m-%d-%Y-%H-%M-%S")))
