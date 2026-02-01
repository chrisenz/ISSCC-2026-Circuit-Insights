# Redefine the log functions
import numpy as np
from numpy import log as ln
from numpy import log10 as log
from numpy import sqrt as sqrt
from numpy import exp as exp
from numpy import tanh as tanh
from numpy import arctan as atan

# Interpolation functions
def vdinter1(vd,vdsat,epsilon):
    return (vd+vdsat-sqrt((vd-vdsat)**2+4*epsilon**2))/2

def vdinter(vd,vdsat,epsilon):
    return vdinter1(vd,vdsat,epsilon)-vdinter1(0,vdsat,epsilon)

def vdss0(vd,vdsat,epsilon):
    return (vd-vdsat+sqrt((vd-vdsat)**2+4*epsilon**2))/2

def vdss(vd,vdsat,epsilon):
    return vdss0(vd,vdsat,epsilon)-vdss0(0,vdsat,epsilon)

# Velocity-field models
def velsat1(e):
    if e<=1:
        vdrift=e
    else:
        vdrift=1
    return vdrift

def velsat2(e,a):
    return e/pow(1+pow(e,a),1/a)

# Basic EKV functions
def vp_vg(vg,vt0,n):
    return (vg-vt0)/n

from scipy.special import lambertw
def q_v_lambert(v):
    return np.real(lambertw(2*exp(v))/2)

def q_v(v):
    if v <= -15:
        q0=1/(2+exp(-v))
    else:
        if v <-0.35:
            z0=1.55+exp(-v)
        else:
            z0=2/(1.3+v-ln(v+1.6))
        z1=(2+z0)/(1+v+ln(z0))
        q0=(1+v+ln(z1))/(2+z1)
    return q0

def q_vp_vx(vpx,vt0,n):
    vp=vp_vg(vg,vt0,n)
    return q_v(vp-vx)

def i_q(q):
    return q**2+q

def id_q(qs,qd):
    ifor=i_q(qs)
    irev=i_q(qd)
    return ifor-irev

def id0_v(vg,vs,vd,vt0,n):
    vp=vp_vg(vg,vt0,n)
    qs=q_v(vp-vs)
    qd=q_v(vp-vd)
    return id_q(qs,qd)

# General velocity saturation functions
def vdsat_v_q(vp,qdsat):
    return vp-2*qdsat-ln(qdsat)

def idsat1_qs(qs,lc):
    return 2*qdsat1_qs(qs,lc)/lc

def idsat_qdsat(qdsat,lc):
    return 2*qdsat/lc

# Velocity saturation model 1
def qdsat1_qs(qs,lc):
    return 2*lc*(qs**2+qs)/(2+lc+sqrt(4*(1+lc)+lc**2*(1+2*qs)**2))

def vdsat1_vg_vs(vg,vs,vt0,n,lc):
    vp=vp_vg(vg,vt0,n)
    qs=q_vg_vx(vg,vs,vt0,n)
    qdsat=qdsat1_qs(qs,lc)
    return vdsat(vp,qdsat)

def qd_prime1_v(vg,vs,vd,vt0,n,lc):
    qs=q_vg_vx(vg,vs,vt0,n)
    vdsat=vdsat_vg_vs(vg,vs,vt0,n,lc)
    if vd<=vdsat:
        qdprime=q_vg_vx(vg,vd,vt0,n)
    else:
        qs=q_vg_vx(vg,vs,vt0,n)
        qdprime=qdsat1_qs(qs,lc)
    return qdprime

def qd_prime1_q(qd,qdsat):
    if qd>qdsat:
        qdprime=qd
    else:
        qdprime=qdsat
    return qdprime

def id1_v(vg,vs,vd,vt0,n,lc):
    vp=vp_vg(vg,vt0,n)
    qs=q_v(vp-vs)
    qd=q_v(vp-vd)
    qdsat=qdsat1_qs(qs,lc)
    qdprime=qd_prime1_q(qd,qdsat)
    return id_q(qs,qdprime)

def id2_v(vg,vs,vd,vt0,n,lc,gds):
    vp=vp_vg(vg,vt0,n)
    qs=q_v(vp-vs)
    qd=q_v(vp-vd)
    qdsat=qdsat1_qs(qs,lc)
    vdsat=vdsat_v_q(vp,qdsat)
    idsat=2*qdsat/lc
    if vd<=vdsat:
        qdprime=qd_prime1_q(qd,qdsat)
        id=id_q(qs,qdprime)
    else:
        id=idsat+gds*(vd-vdsat)
    return id

def id3_v(vd,vdsat,gds):
    if vd > vdsat:
        id=gds*(vd-vdsat)
    else:
        id=0
    return id

def id4_v(vg,vs,vd,vt0,n,lc,gds,vdsat):
    vp=vp_vg(vg,vt0,n)
    qs=q_v(vp-vs)
    qd=q_v(vp-vd)
    qdsat=qdsat1_qs(qs,lc)
    idsat=2*qdsat/lc
    if vd<=vdsat:
        id=id_q(qs,qdprime)
    else:
        id=idsat+gds*(vd-vdsat)
    return id

# Velocity saturation model 2
def id_q_mod2(qs,qd,lc):
    id0=id_q(qs,qd)
    return id0/(1+lc*(qs-qd))

def id_v_mod2(vg,vs,vd,vt0,n,lc):
    vp=vp_vg(vg,vt0,n)
    qs=q_v(vp-vs)
#    qs=q_v_lambert(vp-vs)
    qd=q_v(vp-vd)
#    qd=q_v_lambert(vp-vd)
    return id_q_mod2(qs,qd,lc)

def qdpeak_qs(qs,lc):
    return (1+lc*qs-sqrt(1+lc*(1+2*qs)))/lc

#def qdsat2_qs(qs,lc):
#    return 2*lc*(qs**2+qs)/(2+lc*(1+2*qs)+sqrt(4+lc**2+4*lc*(1+2*qs)**2))

def qdsat2_qs(qs,lc):
    return (2+lc*(1+2*qs)-sqrt(4+lc*(4+lc+8*qs))/(2*lc))

def qd_prime2_v(vp,vd,vdsat,qdsat):
    if vd<=vdsat:
        qdprime=q_v(vp-vd)
    else:
        qdprime=qdsat
    return qdprime

def gmd2_q(qs,qd,lc):
    return qd*(1+2*qd-lc*(qs-qd)**2)/((1+2*qd)*(1-lc*(qs-qd))**2)

def qdlim_q(qs,lc,gds):
    return ((1+2*lc*gds)*(1+lc*qs)-sqrt((1+2*lc*gds)*(1+2*lc*qs)))/(lc*(1+2*lc*gds))

def qdlim2_q(qs,lc,gds):
    return ((1+2*lc*gds)*(1+lc*qs)-sqrt((1+2*lc*gds)*(1+lc*(1+2*qs))))/(lc*(1+2*lc*gds))

