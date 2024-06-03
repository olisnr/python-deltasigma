from __future__ import division
from ._LCparam2tf import LCparam2tf

def LCoptparam2tf(x, param):
    """Convert optimization parameters to transfer functions."""
    n = param['n']
    param.update({'gw':np.ones(n-1)})
    if param['form'] == 'FB':
        param.update({'gu':np.concatenate((np.array([1]), np.zeros(n)))})
        param.update({'gv':copy.copy(x)})
        param.update({'gx':np.concatenate((np.zeros(n-1), np.array([1])))
        param.update({'rx':np.zeros(n)})
    elif param['form'] == 'FF':
        param.update({'gu':np.concatenate((np.array([1]), np.zeros(n-1),
                                           np.array([1])))})
        param.update({'gv':np.concatenate((np.array([1]), np.zeros(n))})
        param.update({'gx':copy.copy(x)
        param.update({'rx':np.zeros(n)})
    elif param['form'] == 'FFR':
        param.update({'gu':np.concatenate((np.array([1]), np.zeros(n-1),
                                           np.array([1])))})
        param.update({'gv':np.concatenate((np.array(x[0]), np.zeros(n)))})
        param.update({'gx':np.concatenate((np.zeros(n-1), np.array([1])))})
        param.update({'rx':np.concatenate((np.array([0]), x[1:n]))})
    elif param['form'] == 'GEN':
        param.update({'gv':np.atleast_1d(x[:n])})
        param.update({'rx':np.concatenate((np.array([0]), x[n:2*n-1]))})
    else:
        raise ValueError('form=%s is not supported.' % form)
    H, L0, ABCD = LCparam2tf(param)
    return H, L0, ABCD, param
