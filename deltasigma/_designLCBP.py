from __future__ import division

import numpy as np
import copy

def designLCBP6(n=3, OSR=64, opt=2, Hinf=1.6, f0=0.25, t=None, form='FB', x0=None, dbg=False):
    param = {'n':n, 'OSR':OSR, 'Hinf':Hinf, 'f0':f0, 't':t, 'form':form}
    if opt:
        w = 2*np.pi*(f0 + 0.25/OSR*ds_optzeros(n, opt))
    else:
        w = 2*np.pi*f0*np.ones(n)
    param.update({'l':1.0/w})
    param.update({'c':1.0/w})
    if x0 is None:
        x = np.arange(2, n+1)
    else:
        if len(x0) != n:
            raise ValueError('Initial guess (x0) has the wrong number of parameters.')
        else:
            x = copy(x0)
    if dbg:
        fprintf(1,char('Iterations for initial pole placement:\n'))
    options=optimset(char('TolX'),0.001,char('TolFun'),1,char('MaxIter'),1000)
    options=optimset(options,char('LargeScale'),char('off'))
    x=fmincon(LCObj1a,x,[],[],[],[],[],[],LCObj1b,options,param,0.97,dbg)
    H=LCoptparam2tf(x,param)
    rmax=max(abs_(H.p[:]))
    if rmax > 0.97:
        fprintf(2,char('Unable to find parameter values which stabilize the system.\n'))
        return param,H,L0,ABCD,x
    if dbg:
        fprintf(1,char('\nParameter Optimization:\n'))
    options=optimset(char('TolX'),0.0001,char('TolFun'),0.5,char('TolCon'),0.01)
    options=optimset(options,char('MaxIter'),1000)
    options=optimset(options,char('DiffMaxChange'),0.01)
    options=optimset(options,char('DiffMinChange'),0.0001)
    options=optimset(options,char('LargeScale'),char('off'))
    if dbg > 1:
        options=optimset(options,char('Display'),char('iter'))
    x=fmincon(LCObj2a,x,[],[],[],[],[],[],LCObj2b,options,param,0.97,dbg)
    H,L0,ABCD,param=LCoptparam2tf(x,param,nargout=4)
    gain_f0=abs_(evalTFP(L0,H,f0))
    param.gu = copy(param.gu / gain_f0)
    L0.k = copy(L0.k / gain_f0)
    ABCD[:,2 * n + 1]=ABCD[:,2 * n + 1] / gain_f0
    if strcmp(form,char('FF')) or strcmp(form,char('FFR')):
        param.gu[n + 1]=1
    if dbg:
        LCplotTF(H,L0,param)
    return param,H,L0,ABCD,x

def _LCObj1a(x=None, param=None, max_radius=None, dbg=None):
    """The objective function for the initial optimization
    process used to put the roots of the denominator of the LCBP NTF inside
    the unit circle."""
    return 1. # No actual objective

def LCObj1b(x=None,param=None,max_radius=None,dbg=None,*args,**kwargs):
    """The constraints function for the initial optimization
    process used to put the roots of the denominator of the LCBP NTF inside
    the unit circle."""
    H = LCoptparam2tf(x,param)
    objective = 1
    rmax=max(abs_(H.p[:]))
    C = rmax - max_radius
    Ceq = np.array([])
    if dbg:
        print('x = [ %.4f ]\n' % x)
        print('rmax = %f\n\n' % rmax)
    return C, Ceq

def LCObj2a(x=None,param=None,rmax=None,dbg=None,*args,**kwargs):
    nargout = kwargs["nargout"] if kwargs else None
    varargin = cellarray(args)
    nargin = 4-[x,param,rmax,dbg].count(None)+len(args)

    if dbg:
        fprintf(1,char('x = ['))
        fprintf(1,char('%8.5f '),x)
        fprintf(1,char(']\n'))
    H,L0,ABCD=LCoptparam2tf(x,param,nargout=3)
    f1=param.f0 - 0.25 / param.OSR
    f2=param.f0 + 0.25 / param.OSR
    objective=dbv(rmsGain(H,f1,f2) / sqrt(3 * param.OSR))
    fprintf(1,char('N0^2 = %.1fdB\n'),objective)
    return objective

def LCObj2b(x=None,param=None,rmax=None,dbg=None,*args,**kwargs):
    nargout = kwargs["nargout"] if kwargs else None
    varargin = cellarray(args)
    nargin = 4-[x,param,rmax,dbg].count(None)+len(args)

    H,L0,ABCD=LCoptparam2tf(x,param,nargout=3)
    Ceq=matlabarray([])
    max_radius=max(abs_(H.p[:]))
    H_inf=infnorm(H)
    stf0=abs_(evalTFP(L0,H,param.f0))
    tmp1,tmp2,tmp3,y=simulateDSM(0.5 / stf0 * sin(2 * pi * param.f0 * cat(arange(0,1000))),ABCD,nargout=4)
    ymax=max(abs_(y))
    if strcmp(param.form,char('FB')):
        f1=param.f0 - 0.25 / param.OSR
        f2=param.f0 + 0.25 / param.OSR
        stf1=abs_(evalTFP(L0,H,f1)) / stf0
        stf2=abs_(evalTFP(L0,H,f2)) / stf0
        C=matlabarray(cat(100 * (max_radius - rmax),H_inf - param.Hinf,stf1 - 1.01,stf2 - 1.01(ymax - 5) / 10))
    else:
        C=matlabarray(cat(100 * (max_radius - rmax),H_inf - param.Hinf(ymax - 5) / 10))
    if dbg:
        fprintf(1,char('constraints ='))
        fprintf(1,char(' %8.5f'),C)
        fprintf(1,char('\n'))
        fprintf(1,char('rmax=%5.3f, Hinf = %4.2f, ymax = %5.3f\n\n'),max_radius,H_inf,ymax)
    return C,Ceq
