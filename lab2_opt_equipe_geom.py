import time

import numpy as np
from scipy.optimize import minimize

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

rad2deg = 180/np.pi

DV = [('AR_w',      8.0,   7.0,  12.0),
      ('xr_w',     17.0,  14.0,  20.0),
      ('S_w',     390.0, 330.0, 460.0),
      ('sweep_w',  0.58,  0.40,  0.70),
      ('x_mlg',    31.0,  28.0,  35.0),
      ('tcr_w',     0.18,  0.12,  0.24)]

dv_names = [d[0] for d in DV]
Xref = np.array([d[1] for d in DV])
bounds_phys = np.array([[d[2], d[3]] for d in DV])

CON = [('deltaS_wlan',      '>=',  0.0),
       ('SM_fwd',           '<=',  0.30),
       ('SM_aft',           '>=',  0.05),
       ('frac_nlg_fwd',     '<=',  0.18),
       ('frac_nlg_aft',     '>=',  0.03),
       ('alpha_tipback',    '>=', 15.0),
       ('alpha_tailstrike', '>=', 10.0),
       ('phi_overturn',     '<=', 63.0),
       ('tank_excess',      '>=',  0.0),
       ('b_w',              '<=', 65.0),
       ('mlg_track',        '<=', 14.0),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

con_names = [c[0] for c in CON]

Xlist = []
flist = []
glist = []

neval = [0]

cache = {'x': None, 'out': None}

def run_analysis(x):

    if cache['x'] is not None and np.array_equal(x, cache['x']):
        return cache['out']

    X = np.asarray(x)*Xref

    airplane = standard_airplane('my_airplane')
    for name, value in zip(dv_names, X):
        airplane['inputs'][name] = value

    analyze(airplane)
    neval[0] = neval[0] + 1

    g = airplane['geometry']
    b = airplane['balance']
    lg = airplane['landing_gear']
    tm = airplane['thrust_matching']
    i = airplane['inputs']

    # corda e bordo de fuga da asa na estacao lateral do trem principal
    eta = i['y_mlg']/(g['b_w']/2)
    c_mlg = g['cr_w'] - (g['cr_w'] - g['ct_w'])*eta
    te_mlg = i['xr_w'] + i['y_mlg']*np.tan(i['sweep_w']) + c_mlg

    out = {'W0'               : tm['W0']/gravity,
           'Wf'               : tm['W_fuel']/gravity,
           'deltaS_wlan'      : tm['deltaS_wlan']/i['S_w'],
           'SM_fwd'           : b['SM_fwd'],
           'SM_aft'           : b['SM_aft'],
           'frac_nlg_fwd'     : lg['frac_nlg_fwd'],
           'frac_nlg_aft'     : lg['frac_nlg_aft'],
           'alpha_tipback'    : lg['alpha_tipback']*rad2deg,
           'alpha_tailstrike' : lg['alpha_tailstrike']*rad2deg,
           'phi_overturn'     : lg['phi_overturn']*rad2deg,
           'tank_excess'      : b['tank_excess'],
           'b_w'              : g['b_w'],
           'mlg_track'        : lg['mlg_track'],
           'h_tail'           : (i['zr_v'] + g['b_v']) - i['z_lg'],
           'T0req'            : max(tm['T0req'].values())/tm['T0'],
           # raiz da empenagem vertical deve terminar dentro da fuselagem
           'vt_fit'           : (g['xr_v'] + g['cr_v'])/i['L_f'],
           # trem principal deve ficar sob a asa, nao atras dela
           'mlg_fit'          : i['x_mlg']/te_mlg}

    cache['x'] = np.array(x)
    cache['out'] = out

    Xlist.append(X)
    flist.append(out['W0'])
    glist.append(constraints(x))

    return out

def objfun(x):

    out = run_analysis(x)

    f = out['W0']/W0_ref

    return f

def constraints(x):

    out = run_analysis(x)

    g = []
    for name, sense, lim in CON:
        v = out[name]
        if lim == 0.0:
            g.append(v if sense == '>=' else -v)
        elif sense == '>=':
            g.append(v/lim - 1)
        else:
            g.append(1 - v/lim)

    return np.array(g)

cons = [{'type': 'ineq', 'fun': constraints}]

X0 = Xref.copy()
x0 = X0/Xref
bounds = np.sort(bounds_phys/Xref[:,None], axis=1)

out_ref = run_analysis(x0)
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']
g_ref = constraints(x0)

Xlist.clear()
flist.clear()
glist.clear()
neval[0] = 0
cache['x'] = None

options = {'maxiter': 200, 'ftol': 1e-8}

t_start = time.time()
result = minimize(objfun, x0,
                  constraints = cons, bounds = bounds,
                  method = 'slsqp', options = options)
t_elapsed = time.time() - t_start

print(result)

xopt = result.x
Xopt = xopt*Xref
out_opt = run_analysis(xopt)
g_opt = constraints(xopt)

print('')
print('='*70)
print('VARIAVEIS DE PROJETO')
print('='*70)
print('%-18s %12s %12s %12s'%('', 'inicial', 'otimizado', 'variacao'))
for k, name in enumerate(dv_names):
    print('%-18s %12.4f %12.4f %11.2f%%'%(name, X0[k], Xopt[k],
                                          100*(Xopt[k]-X0[k])/X0[k]))

print('')
print('='*70)
print('OBJETIVO')
print('='*70)
print('%-18s %12.1f %12.1f %11.2f%%'%('W0 [kgf]', W0_ref, out_opt['W0'],
                                      100*(out_opt['W0']-W0_ref)/W0_ref))
print('%-18s %12.1f %12.1f %11.2f%%'%('Wf [kgf]', Wf_ref, out_opt['Wf'],
                                      100*(out_opt['Wf']-Wf_ref)/Wf_ref))

print('')
print('='*70)
print('RESTRICOES  (g normalizada >= 0; ATIVA se |g| < 1e-4)')
print('='*70)
print('%-18s %10s %10s %10s %10s  %s'%('', 'inicial', 'otimizado', 'limite', 'g_norm', 'estado'))

tol_ativo = 1e-4
n_ativas = 0
for k, (name, sense, lim) in enumerate(CON):
    if abs(g_opt[k]) < tol_ativo:
        estado = 'ATIVA'
        n_ativas += 1
    elif g_opt[k] < 0:
        estado = 'VIOLADA'
    else:
        estado = ''
    print('%-18s %10.4f %10.4f %4s %5.2f %10.2e  %s'%(name, out_ref[name], out_opt[name],
                                                      sense, lim, g_opt[k], estado))

print('')
for k, name in enumerate(dv_names):
    for b in bounds[k]:
        if abs(xopt[k]-b) < tol_ativo:
            print('bound ATIVO: %s = %.4f'%(name, b*Xref[k]))
            n_ativas += 1

print('')
print('='*70)
print('DESEMPENHO')
print('='*70)
print('Iteracoes do SLSQP (nit):      %d'%result.nit)
print('Chamadas de objfun (nfev):     %d'%result.nfev)
print('Execucoes do analyze:          %d'%neval[0])
print('Tempo de otimizacao:           %.2f s'%t_elapsed)
print('Restricoes ativas:             %d'%n_ativas)
print('='*70)
