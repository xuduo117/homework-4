#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 13:12:44 2017

@author: xuduo
"""


import numpy as np
import matplotlib.pyplot as plt
import corner
import matplotlib
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import scipy
    
def fit_function(a_0,a_1,a_2,x):
    return -(2*np.pi*G)**(1./3.)*(1.1*M_sun)**(-2./3.)*a_0*M_jup*(a_1*day)**(-1./3.)*np.sin((x-a_2)*2*np.pi/(a_1))/1e2

def phase_function(a_0,a_1,a_2,x):
    return -(2*np.pi*G)**(1./3.)*(1.1*M_sun)**(-2./3.)*a_0*M_jup*(a_1*day)**(-1./3.)*np.sin(((x)*2*np.pi/(a_1))%(2*np.pi))/1e2


def log_gauss_function(x1,x2,sigma):
    return -(x1-x2)**2/2.0/sigma**2


def log_likelihood(a0,a1,a2):
    return np.sum(log_gauss_function(data_2_y,fit_function(a0,a1,a2,data_2_x),data_2_sigma))

def alpha(a0_c,a1_c,a2_c,a0,a1,a2):
    return np.min([1.0,np.exp(log_likelihood(a0_c,a1_c,a2_c)-log_likelihood(a0,a1,a2))])


def log_likelihood_array_a0(a0,a1,a2):
    llh_a0=np.array([])
    for ctt_1 in range(len(a0)):
        llh=np.sum(log_gauss_function(data_2_y,fit_function(a0[ctt_1],a1,a2,data_2_x),data_2_sigma))
        llh_a0=np.append(llh_a0,llh)
    return llh_a0

def log_likelihood_array_a1(a0,a1,a2):
    llh_a1=np.array([])
    for ctt_1 in range(len(a1)):
        llh=np.sum(log_gauss_function(data_2_y,fit_function(a0,a1[ctt_1],a2,data_2_x),data_2_sigma))
        llh_a1=np.append(llh_a1,llh)
    return llh_a1

def log_likelihood_array_a2(a0,a1,a2):
    llh_a2=np.array([])
    for ctt_1 in range(len(a2)):
        llh=np.sum(log_gauss_function(data_2_y,fit_function(a0,a1,a2[ctt_1],data_2_x),data_2_sigma))
        llh_a2=np.append(llh_a2,llh)
    return llh_a2

#def function_cnst(x,max_x):
#    return max_x
# 
#def int_function_cnst(x,max_x):
#    return max_x*x
#
#def invint_function_cnst(x,max_x):
#    return 1./max_x*x
#    
#def next_step(a2,llh_a2):
#    
#    prob_a2=np.exp(llh_a2-np.max(llh_a2))
#    indx_a2=np.argwhere(prob_a2>0)
##    f=interp1d(a2,prob_a2,kind='linear',fill_value=0)
#
#    xmin=a2[indx_a2[0]]
#    xmax=a2[indx_a2[-1]]
#    
#    ran=np.random.uniform(low=xmin,high=xmax,size=2)
#    
##    n=0
##    max_prob=np.max(prob_a2)*1.1
##    
##    #"""
##    x=np.linspace(xmin,xmax,1000)  
##    y=f(x)  
##    pmin=0.  
##    pmax=int_function_cnst(xmax,max_prob)
##    #pmax=1
##       
##    naccept=0  
##    ntrial=0  
##       
##    ran=[] 
##    while naccept==n:  
##        y=np.random.uniform(pmin,pmax) 
##        x0=invint_function_cnst(y,max_prob)
##        f1=np.random.uniform(0,function_cnst(x0,max_prob))
##        
##        ntrial=ntrial+1 
##        if f1<f(x0):
##            ran.append(x0)  
##            naccept=naccept+1  
#    return ran[0]


def plot_hist1d_1(a0_all,title,scinot):
    a0_median=np.median(a0_all)
    a0_16=np.percentile(a0_all,16)
    a0_84=np.percentile(a0_all,84)
    print a0_median,a0_16-a0_median,a0_84-a0_median
    plt.clf()
    n=plt.hist(a0_all,bins=30)
    plt.plot([a0_16,a0_16], [0,np.max(n[0])*1.1],linestyle='dashed')
    plt.plot([a0_median,a0_median], [0,np.max(n[0])*1.1],linestyle='dashed')
    plt.plot([a0_84,a0_84], [0,np.max(n[0])*1.1],linestyle='dashed')
    plt.title(title)
    plt.xlabel(title)
    if scinot !=0:
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.text(a0_16,np.max(n[0])*0.8,str('{:.5f}'.format(a0_16)),horizontalalignment='center',fontsize=13)
        plt.text(a0_median,np.max(n[0])*0.9,str('{:.5f}'.format(a0_median)),horizontalalignment='center',fontsize=13)
        plt.text(a0_84,np.max(n[0])*1.0,str('{:.5f}'.format(a0_84)),horizontalalignment='center',fontsize=13)
    else:
        plt.text(a0_16,np.max(n[0])*0.8,str('{:.4f}'.format(a0_16)),horizontalalignment='center',fontsize=13)
        plt.text(a0_median,np.max(n[0])*0.9,str('{:.4f}'.format(a0_median)),horizontalalignment='center',fontsize=13)
        plt.text(a0_84,np.max(n[0])*1.0,str('{:.4f}'.format(a0_84)),horizontalalignment='center',fontsize=13)
    plt.ylabel('N')
    plt.savefig(str(title)+'_1d.pdf',bbox_inches='tight')


#def find_confi_level(llh_all_sort,level):
#    array_len=len(llh_all_sort)
#    num=array_len/2
#    n=1
#    while np.sum(llh_all_sort[0:num])>level or np.sum(llh_all_sort[0:num+1])<level:
#        if np.sum(llh_all_sort[0:num+1])<level:
#            num_new=num+np.floor(num/2.0/n)
#        if np.sum(llh_all_sort[0:num])>level:
#            num_new=num-np.ceil(num/2.0/n)
#        num=np.int(num_new)
#        n=n+1
#    return num 

def find_confi_level_2(llh_all_sort,level):
    array_len=len(llh_all_sort)
    num=array_len/2
    n=1
    while np.sum(llh_all_sort[0:num])>level or np.sum(llh_all_sort[0:num+1])<level:
        if np.sum(llh_all_sort[0:num+1])<level:
            num_new=num+np.floor(num/2.0/n)
        if np.sum(llh_all_sort[0:num])>level:
            num_new=num-np.ceil(num/2.0/n)
        num=np.int(num_new)
        n=n+1
    return num 


#def plot_hist2d_1(a0_all,a1_all,xtitle,ytitle):
#    x_grid = np.linspace(min(a0_all),max(a0_all),100)
#    y_grid = np.linspace(min(a1_all),max(a1_all),100)
#            
#    z_grid = griddata((a0_all, a1_all), llh_all/np.sum(llh_all), (x_grid[None,:], y_grid[:,None]), method='nearest')
#    
#    z_grid=scipy.ndimage.gaussian_filter(z_grid,3, output=None,mode='reflect', cval=0.0, truncate=4.0)
#    plt.clf()
#    
#    plt.hist2d(a0_all, a1_all,bins=40)
#    cb=plt.colorbar()
#    cb.ax.tick_params(labelsize=14)
#    #levels=[llh_all_68,llh_all_90,llh_all_99]
#    levels=[llh_all_99,llh_all_90,llh_all_68]
#    CS=plt.contour(x_grid,y_grid,z_grid,levels=levels,linewidths=1.5, colors='white')
#    
#    fmt={}
#    #strs=['68%','90%','99%']
#    strs=['99.7%','95%','68%']
#    for l, s in zip(CS.levels, strs):
#        fmt[l] = s
#    
#    plt.clabel(CS,CS.levels[:], inline=True,fmt=fmt, fontsize=14,thick=20)
#    plt.xlabel(xtitle,fontsize=14)
#    plt.ylabel(ytitle,fontsize=14)
#    plt.savefig(str(xtitle)+'_'+str(ytitle)+'_2d.pdf',bbox_inches='tight')
#    plt.clf()


def plot_hist2d_2(a0_all,a1_all,xtitle,ytitle):
    plt.clf()

    counts, xedges, yedges, Image=plt.hist2d(a0_all, a1_all,bins=40)
    xedges_1=xedges[0:-1]+(xedges[1]-xedges[0])/2.0
    yedges_1=yedges[0:-1]+(yedges[1]-yedges[0])/2.0
    counts=scipy.ndimage.gaussian_filter(counts,1.5, output=None,mode='reflect', cval=0.0, truncate=4.0)
    
    cb=plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    #levels=[llh_all_68,llh_all_90,llh_all_99]
    counts_sort=np.sort(counts.flatten())
    counts_sort=counts_sort[::-1]
    counts_sort_norm=counts_sort/np.sum(counts_sort)
    
    counts_68=counts_sort[find_confi_level_2(counts_sort_norm,0.68)]
    counts_95=counts_sort[find_confi_level_2(counts_sort_norm,0.95)]
    counts_99=counts_sort[find_confi_level_2(counts_sort_norm,0.997)]
    
    levels=[counts_99,counts_95,counts_68]
    CS=plt.contour(xedges_1,yedges_1,counts.transpose(),levels=levels,linewidths=1.5, colors='white')
    
    fmt={}
    #strs=['68%','90%','99%']
    strs=['99.7%','95%','68%']
    for l, s in zip(CS.levels, strs):
        fmt[l] = s
    
    plt.clabel(CS,CS.levels[:], inline=True,fmt=fmt, fontsize=14,thick=20)
    plt.xlabel(xtitle,fontsize=14)
    plt.ylabel(ytitle,fontsize=14)
    plt.savefig(str(xtitle)+'_'+str(ytitle)+'_2dhist.pdf',bbox_inches='tight')
    plt.clf()

#"""

G=6.674079999999999e-08
R_sun=69570000000.0
R_jup=7149200000.0
M_sun=1.9884754153381438e+33
M_earth=5.972364730419773e+27
M_jup=1.8981871658715508e+30
AU=14959787070000.0
yr=3.15576e7
day=86400.0



data=np.loadtxt('../HD209458_3_KECK.vels')
xx=data[:,0]-2451341
yy=data[:,1]
y_err=data[:,2]

data_2_x=xx
data_2_y=yy
data_2_sigma=y_err
#plt.errorbar(xx,yy,yerr=y_err)



#guess = np.array([0.685,3.524,-0.046615000103557236])
#guess = np.array([0.65,3.5239,1.7155499998964387])
guess = np.array([0.9,3.524,1.6])

#"""
a0=guess[0]
a1=guess[1]
a2=guess[2]

#print likelihood_1(a0,a1)
accept_num=0
rej_num=0
total_num=0
a0_all=[]
a1_all=[]
a2_all=[]

#plt.figure(0)
scale_normal=1.0  ##  0.04  1.0







#plt.hist(np.asarray(ran),normed=True)
#plt.plot()




while accept_num < 100000 and total_num <400000:
    a0_c=np.abs(np.random.normal(loc=a0,scale=scale_normal*1e-3))
    a1_c=np.abs(np.random.normal(loc=a1,scale=scale_normal*4e-5))
    a2_c=np.random.normal(loc=a2,scale=scale_normal*4e-4)

    
#    a0_array=np.linspace(1e-4,2,200)
#    llh_a0=log_likelihood_array_a0(a0_array,a1,a2)
#    
#    a0_c=next_step(a0_array,llh_a0)
#    
#    a1_array=np.linspace(1e-4,10,200)
#    llh_a1=log_likelihood_array_a1(a0_c,a1_array,a2)
#    
#    a1_c=next_step(a1_array,llh_a1)
#    
#    a2_array=np.linspace(1e-4,3.5,200)
#    llh_a2=log_likelihood_array_a2(a0_c,a1_c,a2_array)
#    
#    a2_c=next_step(a2_array,llh_a2)
    
    u=np.random.uniform(low=0,high=1.0)
    total_num=total_num+1
    if u < alpha(a0_c,a1_c,a2_c,a0,a1,a2):
#    if log_likelihood(a0,a1) < log_likelihood(a0_c,a1_c):            
        a0=a0_c
        a1=a1_c
        a2=a2_c
        accept_num=accept_num+1
    else:
        a0=a0
        a1=a1
        a2=a2
        rej_num=rej_num+1
#    while accept_num in range(0,100,10):
#        plt.plot([a0],[a1],'o')
    a0_all.append(a0)
    a1_all.append(a1)
    a2_all.append(a2)
    
#print a0,a1

#plt.plot(a0_all,a1_all)
#plt.plot(a0_all[::50],a1_all[::50])

print accept_num/(rej_num+accept_num+0.0),accept_num,rej_num+accept_num
samples=np.array([a0_all,a1_all,a2_all]).transpose()
#m_true = 1.5
#b_true = -60

samples=samples[2000:]

a0_median=np.median(samples[:,0])
a1_median=np.median(samples[:,1])
a2_median=np.median(samples[:,2])
#np.save('prob5.1_sample.npy',samples)
#samples=np.load('prob5.1_sample.npy')

matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['xtick.labelsize'] = 13 
matplotlib.rcParams['ytick.labelsize'] = 13 

plt.figure(1)
#plt.clf()x
fig = corner.corner(samples, labels=["$a_0$", "$a_1$","$a_2$"],quantiles=[0.16, 0.5, 0.84],title_fmt=".6f",
                       show_titles=True,levels=(0.68,0.95,0.997),
                       title_kwargs={"fontsize": 16})
#plt.show()
plt.savefig('corner_1.pdf',bbox_inches='tight')


plt.figure(2)
plt.clf()
gs1 = matplotlib.gridspec.GridSpec(3, 3)
#gs1.update(left=0.05, right=0.48, wspace=0.05)
ax1 = plt.subplot(gs1[:-1, :])
ax2 = plt.subplot(gs1[-1, :],sharex = ax1)

ax1.errorbar(xx,yy,yerr=y_err,linestyle='none',marker='o',label='observed')
ax1.plot(xx,fit_function(a0_median,a1_median,a2_median,xx),label='fit',linestyle='dashed')
plt.setp(ax1.get_xticklabels(), visible=False)
ax2.errorbar(xx, yy-fit_function(a0_median,a1_median,a2_median,xx),yerr=y_err)
ax1.set_ylabel('RV (m/s)')
ax2.set_ylabel('residual (m/s)')
ax2.set_xlabel('Time (JD +2451341)')
ax1.legend()
plt.subplots_adjust(hspace=.0)
plt.savefig('fitplot.pdf',bbox_inches='tight')

#plt.clf()
#xx_plot=xx##np.linspace(0,210,num=10000)
#plt.errorbar(xx,yy,yerr=y_err,linestyle='none',marker='o',label='observed')
#plt.plot(xx_plot,fit_function(a0,a1,a2,xx_plot),label='fit',linestyle='--')
#plt.legend()

#"""

plt.figure(3)
plt.clf()
plt.errorbar(((xx-a2_median)%a1_median)/a1_median ,yy,yerr=y_err,linestyle='None',marker='o',label='observed')
x_plot=np.linspace(0,a1_median,1000)
plt.plot(x_plot/a1_median,phase_function(a0_median,a1_median,a2_median,x_plot),label='fit')
plt.xlabel('phase')
plt.ylabel('RV (m/s)')
plt.legend()
plt.savefig('RV_phase.pdf',bbox_inches='tight')


#"""
a0_all=samples[:,0]
a1_all=samples[:,1]
a2_all=samples[:,2]

plot_hist1d_1(a0_all,'m sini',0)
plot_hist1d_1(a1_all,'T',1)
plot_hist1d_1(a2_all,'T0',0)


log_llh_all=np.array([])
for ctt in range(len(a0_all)):
    log_llh_all=np.append(log_llh_all,log_likelihood(a0_all[ctt],a1_all[ctt],a2_all[ctt]))
 
log_llh_all=log_llh_all-np.max(log_llh_all)

llh_all=np.exp(log_llh_all)
llh_all_sort=np.sort(llh_all)
llh_all_sort=llh_all_sort[::-1]
llh_all_sort=llh_all_sort/np.sum(llh_all_sort)




#llh_all_68=llh_all_sort[find_confi_level(llh_all_sort,0.68)]
#llh_all_90=llh_all_sort[find_confi_level(llh_all_sort,0.95)]
#llh_all_99=llh_all_sort[find_confi_level(llh_all_sort,0.997)]

#plot_hist2d_1(a0_all,a1_all,'m sini','T')
#plot_hist2d_1(a0_all,a2_all,'m sini','T0')
#plot_hist2d_1(a1_all,a2_all,'T','T0')


plot_hist2d_2(a0_all,a1_all,'m sini','T')
plot_hist2d_2(a0_all,a2_all,'m sini','T0')
plot_hist2d_2(a1_all,a2_all,'T','T0')





