# This includes all the functions used for the EKV MOSFET model
# Redefine the log functions
import numpy as np
from numpy import log as ln
from numpy import log10 as log
from numpy import sqrt as sqrt
from numpy import exp as exp
from numpy import arctan as atan

#--------------------------------------------------
# Large-signal normalized functions
#--------------------------------------------------

# Normalized drain current versus charge
def iq(q):
    return q*(q+1)

# Normalized charge versus current
def qi(i):
    return (sqrt(4*i+1)-1)/2

# Normalized charge versus inversion coefficient including VS
# New sEKV model
def qs_ic(ic,lc):
    return (sqrt(4*ic+1+(lc*ic)**2)-1)/2

# Normalized saturation voltage vp-v versus the normalized charge
def vq(q):
    return 2*q+ln(q)

def vi(i):
    q=qi(i)
    return 2*q+ln(q)

# Inverse function giving q(v) corresponding to the inverse function of v(q) using the Lambert W-function
from scipy.special import lambertw
def qv(v):
    return np.real(lambertw(2*exp(v))/2)

# Inverse function giving q(v) corresponding to the inverse function of v(q) using the EKV 2.6 approximation
def qvapprox(v):
    if v <= -15:
        q0=1/(2+exp(v))
    else:
        if v <-0.35:
            z0=1.55+exp(-v)
        else:
            z0=2/(1.3+v-ln(v+1.6))
    z1=(2+z0)/(1+v+ln(z0))
    q0=(1+v+ln(z1))/(2+z1)
    return q0

# Normalized saturation voltage vp-vs versus the normalized charge
def vps_ic(ic):
    return ln(sqrt(4*ic+1)-1)+sqrt(4*ic+1)-1-ln(2)

# Inversion coefficient as a function of the normalized saturation voltage vp-vs
def ic_vps(vps):
    return iq(qv(vps))

# Drain-to-source saturation voltage versus inversion coefficient
# New vdssat function to avoid problems in the inverse function (CE: 9.3.2022)
def vdssat_ic(ic):
    vdssatwi=4
    return 2*sqrt(ic+vdssatwi)

# Inversion coefficient versus drain-to-source saturation voltage
# New inverse vdssat function to avoid having negative ic for vdssat<4
def ic_vdssat(vdssat):
    vdssatwi=4
    if vdssat<vdssatwi:
        ic=0
    else:
        ic=(vdssat/2)**2-vdssatwi
    return ic

# Slope factor versus inversion coefficient
def n_ic(ic,gammab,psi0):
    vp=vps_ic(ic)
    return 1+gammab/(2*sqrt(psi0+vp))

#--------------------------------------------------
# Small-signal normalized functions
#--------------------------------------------------

# Source transconductance versus inversion coefficient (long-channel)
def gms_ic_long(ic):
    return (sqrt(4*ic+1)-1)/2

# Normalized Gm/ID function versus inversion coefficient (long-channel)
def gmsid_ic_long(ic):
    return gms_ic_long(ic)/ic

# Normalized source transconductance versus inversion coefficient including velocity saturation
# New sEKV model
def gms_ic(ic,lc):
    num=sqrt(4*ic+1+(lc*ic)**2)-1
    den=2+lc**2*ic
    return num/den

# Normalized Gm/ID function versus inversion coefficient including velocity saturation
# New sEKV model
def gmsid_ic(ic,lc):
    return gms_ic(ic,lc)/ic

# Inverse normalized Gm/ID function (long-channel)
def ic_gmid(gmsid):
    return ((2/gmsid-1)**2-1)/4

#--------------------------------------------------
# Thermal noise factors
#--------------------------------------------------

# Thermal noise excess factor in saturation (long-channel)
def deltan_ic_long(ic):
    qs=qi(ic)
    return 2/3*(qs+3/4)/(qs+1)

# Thermal noise excess factor including short-channel effects
def delta_ic_vs(ic,thetan):
    return delta_ic(ic)*(1+thetan*ic)

# Theoretical thermal noise excess factor including short-channel effects
def gamman_ic_the(ic,n,lambdac):
    gammawi=n/2
    alphan=n/2*lambdac**2
    return gammawi+alphan*ic

#--------------------------------------------------
# Fano noise suppression factor
#--------------------------------------------------

# Fano suppression factor (long-channel)
def fano_ic_long(ic):
    return 2*deltan_ic_long(ic)*gmsid_ic_long(ic)

# Theoretical Fano suppression factor including short-channel effects
def fano_ic_the(ic,n,lambdac):   
    return 2*gamman_ic_the(ic,n,lambdac)/n*gmsid_ic(ic,lambdac)

def fano_ic_cor(ic,lambdac):
    thetan=lambdac**2/(1+lambdac/2)
    return 2*delta_ic_vs(ic,thetan)*gmsid_ic_vs(ic,lambdac)

# Thermal noise excess factor including short-channel effects
def gamma_ic_vs(ic,gammawi,alpha):
    return gammawi+alpha*ic

# Fano suppression factor including short-channel effects
def fano_ic_vs(ic,lambdac,alpha):
    return 2*delta_ic_vs(ic,alpha)*gmsid_ic_vs(ic,lambdac)

def fano_ic_exp(ic,n,lambdac,gammawi,alpha):
    delta=(gammawi+alpha*ic)/n
    return 2*delta*gmsid_ic_vs(ic,lambdac)
