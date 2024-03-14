import os
from math import pi
from auto import *
import numpy as np

def cache_branch(nm, *args, forward=None, backward=None, use_tmp=False, overwrite=True, **kw):
    
    if len(args)==0 or overwrite==False:
        has_for = os.path.exists('b.'+nm+'_for')
        has_back = os.path.exists('b.'+nm+'_back')
        if has_for and has_back:
            return loadbd(nm+'_for') + loadbd(nm+'_back')
        elif has_for:
            return loadbd(nm+'_for')
        elif has_back:
            return loadbd(nm+'_back')
    
    if forward is None:
        forward=kw;
        
    if backward is None:
        backward=kw;
    
    do_forward = forward or (forward == {} and kw=={})
    if do_forward:
        forward = { **kw, **forward}
        print(forward)
        rfor=run(*args, **forward)
        sv(nm+'_for')
        
    do_backward = backward or (backward == {} and kw=={})
    if do_backward:
        backward = { **kw,**{'DS':'-'}, **backward}
        print(backward)
        rback=run(*args, **backward)
        sv(nm+'_back')
    if do_forward and do_backward:
        return rfor + rback
    elif do_forward:
        return rfor
    elif do_backward:
        return rback

        
    

    
    
    
    '''
    Runs some slices for through parameter space for values of a parameter specified in the slices dict()
    
    slices={`par`:[v1, v2, v3,...]}
    
    First we continue in `par`, setting UZ points at all the values v1, v2, v3,...
    
    Then we continue from those UZ points according to ICP.
    
    '''  
def cache_slices(nm, *args, slices={}, overwrite=True, forward=None, backward=None,**kw):
    k=list(slices.keys())[0]
    desired_vals=np.sort(np.array(slices[k]))
    if 'ICP' in kw:
        ICP0=kw['ICP'].copy()
        ICP0[0]=k
    else:
        ICP0=[k]
        
    if 'PAR' in kw:
        k0=kw['PAR'][k]
    else:
        raise KeyError
        
    if 'UZSTOP' not in kw:
        kw['UZSTOP']={}
    
    if np.any(k0>desired_vals) and np.any(k0<desired_vals): #this all assumes ds>0, TODO: negative DS
        UZSTOP0=[desired_vals[0], desired_vals[-1]]
        UZ0 = slices[k][1:-1]
        f0=None
        b0=None
    elif np.any(k0>desired_vals):
        UZSTOP0=[desired_vals[0]]
        UZ0 = slices[k][1:]
        f0=False
        b0=None
    else:
        UZSTOP0=[desired_vals[-1]]
        UZ0 = slices[k][:-1]
        b0=False
        f0=None
        
    kw0={**kw, **{'ICP':ICP0,'UZSTOP':{k:UZSTOP0},'UZR':{k:UZ0}}}
    
    r0 = cache_branch(nm+'_slice0', *args, overwrite=overwrite, forward=f0, backward=b0,  **kw0)
    
    out={v:None for v in desired_vals}
    kw_slices={**kw, **{'PAR':{}}}
    
    actual_vals = [uz.__dict__['b']['data'][0] for uz in r0('UZ') ]
    
    dist = np.abs(np.array(actual_vals)-desired_vals.reshape(-1,1))
    uz_to_slice = {v:i for v,i in zip(actual_vals,np.argmin(dist, axis=0))}
    
    
    for uz, p in zip(r0('UZ'), actual_vals):
        p0 = desired_vals[uz_to_slice[p]]
        out[p0] = cache_branch(f'{nm}_{k}={p0}', uz, ICP=kw['ICP'], UZSTOP=kw['UZSTOP'], overwrite=overwrite, forward=forward, backward=backward )
        
    return out

def bifDiag_measures(bd, m):
    ICP = bd(bd.getLabels()[0]).__dict__['c']['ICP']
    if m not in ICP:
        parnames = bd(bd.getLabels()[0]).__dict__['c']['parnames']
        nms=[p[1] for p in parnames]
        inds=[p[0] for p in parnames]
        i=ICP.index(inds[nms.index(m)])
    else:
        i=ICP.index(m)
        
    if i>=1:
        i=i+1
        
    return np.array([a[i] for a in bd.toArray()])
    

def addLine(p,x,y):
    from auto.graphics import grapher_mpl as grapher 
    
    p.grapher.addArray((x,y))
    grapher.GUIGrapher.plot(p.grapher)
    p.grapher.draw()
