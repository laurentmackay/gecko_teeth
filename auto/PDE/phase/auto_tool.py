import os
from math import pi
from auto import *
import numpy as np

def cache_branch(nm, *args, forward=None, backward=None, use_tmp=False, **kw):
	
    if len(args)==0:
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

		
	

	
	
    
    
def cache_slices(nm, *args, slices={}, **kw):
    '''
    Runs some slices for through parameter space for values of a parameter specified in the slices dict()
    
    slices={`par`:[v1, v2, v3,...]}
    
    First we continue in `par`, setting UZ points at all the values v1, v2, v3,...
    
    Then we continue from those UZ points according to ICP.
    
    '''
    k=list(slices.keys())[0]
    desired_vals=np.sort(np.array(slices[k]))
    ICP0=kw['ICP'].copy()
    ICP0[0]=k
    kw0={**kw, **{'ICP':ICP0},'UZSTOP':{k:[desired_vals[0], desired_vals[-1]]}}

    r0 = cache_branch(nm+'_slice0',*args, UZR=slices[1:-2], **kw0)

    
    out={v:None for v in desired_vals}
    kw_slices={**kw, **{'PAR':{}}}
    

    actual_vals = [uz.__dict__['b']['data'][0] for uz in r0('UZ') ]
    
    uz_to_slice = {v:i for v,i in zip(s,np.argmin(np.abs(np.array(actual_vals)-np.array(t)), axis=0))}

    
    for uz in r0('UZ'):
        p  = uz.__dict__['b']['data'][0]
        p0 = uz_to_slice[p]]
        out[desired_vals[p0]=split_omegas[kappa_close] = cache_branch(f'{nm}_{k}={p0}', uz, ICP=ICP )
        
    return out

    
