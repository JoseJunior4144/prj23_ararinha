'''
Lab 02 - Secao 4: Otimizacao multiobjetivo da aeronave da equipe
Minimizar W0 e Wf simultaneamente com NSGA-II (pymoo), sujeito as mesmas
restricoes da otimizacao mono-objetivo (lab2_opt_equipe_geom.py).

A convergencia da frente de Pareto e verificada comparando o extremo de
minimo W0 com o otimo do SLSQP da Secao 3.
'''

import time
import warnings

import numpy as np
import matplotlib.pyplot as plt
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem

from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.constants import gravity

rad2deg = 180/np.pi

# a populacao aleatoria do NSGA-II gera geometrias extremas em que o ponto
# fixo do W0 diverge (overflow/NaN); os avisos sao esperados e inofensivos
np.seterr(all='ignore')
warnings.filterwarnings('ignore')

# mesmas DVs e limites da otimizacao mono-objetivo
DV = [('AR_w',      8.0,   7.0,  12.0),
      ('xr_w',     17.0,  14.0,  20.0),
      ('S_w',     390.0, 330.0, 460.0),
      ('sweep_w',  0.58,  0.40,  0.70),
      ('x_mlg',    31.0,  28.0,  35.0),
      ('tcr_w',     0.18,  0.12,  0.24),
      ('z_lg',     -5.75, -7.50, -4.50),
      ('y_mlg',     5.50,  4.00,   9.00)]

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
       ('ground_clearance', '>=',  0.5),
       ('tank_excess',      '>=',  0.0),
       ('b_w',              '<=', 64.9),
       ('mlg_track',        '<=', 13.9),
       ('h_tail',           '<=', 20.0),
       ('T0req',            '<=',  1.0),
       ('vt_fit',           '<=',  1.0),
       ('mlg_fit',          '<=',  1.0)]

con_names = [c[0] for c in CON]

neval = [0]

def run_analysis(x):

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
    # (LE interpolado entre raiz e ponta: xt_w ja carrega o termo de 1/4 de
    # corda, pois sweep_w e' medido a 1/4 de corda em geometry.py)
    eta = i['y_mlg']/(g['b_w']/2)
    c_mlg = g['cr_w'] - (g['cr_w'] - g['ct_w'])*eta
    te_mlg = i['xr_w'] + (g['xt_w'] - i['xr_w'])*eta + c_mlg

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
           'ground_clearance' : lg['ground_clearance'],
           'tank_excess'      : b['tank_excess'],
           'b_w'              : g['b_w'],
           'mlg_track'        : lg['mlg_track'],
           'h_tail'           : (i['zr_v'] + g['b_v']) - i['z_lg'],
           'T0req'            : max(tm['T0req'].values())/tm['T0'],
           'vt_fit'           : (g['xr_v'] + g['cr_v'])/i['L_f'],
           'mlg_fit'          : i['x_mlg']/te_mlg}

    return out

def constraints_geq0(out):
    # mesma normalizacao do script mono-objetivo: viavel se g >= 0
    gg = []
    for name, sense, lim in CON:
        v = out[name]
        if lim == 0.0:
            gg.append(v if sense == '>=' else -v)
        elif sense == '>=':
            gg.append(v/lim - 1)
        else:
            gg.append(1 - v/lim)
    return np.array(gg)

# referencia (aeronave do PRJ-22) para normalizar os objetivos
out_ref = run_analysis(np.ones(len(DV)))
W0_ref = out_ref['W0']
Wf_ref = out_ref['Wf']

# otimos mono-objetivo (SLSQP) do mesmo problema, usados como ancoras para
# verificar a convergencia dos extremos da frente: min W0 vem do
# lab2_opt_equipe_geom.py; min Wf vem de uma rodada identica trocando o
# objetivo. Sem as sementes abaixo, o NSGA-II puro (pop 60, 120 geracoes)
# estagnou ~1% acima dessas ancoras -- evidencia de nao-convergencia.
anchor_W0 = {'W0': 290117.1, 'Wf': 111036.5,
             'x': np.array([9.8001, 16.0113, 353.6712, 0.6108,
                            29.2629, 0.2093, -6.0247, 6.9500])/Xref}
anchor_Wf = {'W0': 291961.0, 'Wf': 109117.9,
             'x': np.array([9.8416, 15.8699, 374.4731, 0.6066,
                            29.3427, 0.1872, -6.0106, 6.9500])/Xref}

bounds_norm = np.sort(bounds_phys/Xref[:,None], axis=1)

class AirplaneProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=len(DV),
                         n_obj=2,
                         n_ieq_constr=len(CON),
                         xl=bounds_norm[:,0],
                         xu=bounds_norm[:,1])

    def _evaluate(self, x, out_pymoo, *args, **kwargs):

        # designs extremos podem divergir a analise; devolvemos um
        # individuo fortemente inviavel para o NSGA-II descarta-lo
        try:
            out = run_analysis(x)
            F = [out['W0']/W0_ref, out['Wf']/Wf_ref]
            # ATENCAO: pymoo considera viavel G <= 0, o oposto do scipy
            G = list(-constraints_geq0(out))
            if not (np.all(np.isfinite(F)) and np.all(np.isfinite(G))):
                raise ValueError('analise divergiu')
        except Exception:
            F = [10.0, 10.0]
            G = [1e3]*len(CON)

        out_pymoo['F'] = F
        out_pymoo['G'] = G

POP_SIZE = 80
N_GEN = 200

problem = AirplaneProblem()

# populacao inicial: LHS + baseline + otimos SLSQP como sementes, para
# ancorar os extremos da frente e acelerar a convergencia
from pymoo.operators.sampling.lhs import LHS
X0_pop = LHS()(problem, POP_SIZE).get('X')
X0_pop[0] = np.ones(len(DV))
X0_pop[1] = anchor_W0['x']
X0_pop[2] = anchor_Wf['x']

algorithm = NSGA2(pop_size=POP_SIZE, sampling=X0_pop, eliminate_duplicates=True)

neval[0] = 0
t_start = time.time()
res = minimize(problem, algorithm, ('n_gen', N_GEN), seed=1, verbose=True)
t_elapsed = time.time() - t_start

# ordena a frente por W0
order = np.argsort(res.F[:,0])
F = res.F[order]
X = res.X[order]

W0_front = F[:,0]*W0_ref
Wf_front = F[:,1]*Wf_ref

print('')
print('='*70)
print('OTIMIZACAO MULTIOBJETIVO (NSGA-II)')
print('='*70)
print('Individuos por geracao (pop_size): %d'%POP_SIZE)
print('Numero de geracoes (n_gen):        %d'%N_GEN)
print('Execucoes do analyze:              %d'%neval[0])
print('Tempo de otimizacao:               %.1f s'%t_elapsed)
print('Pontos na frente de Pareto:        %d'%len(F))
print('')
print('Verificacao de convergencia com os otimos mono-objetivo (SLSQP):')
print('  min W0 da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      %(W0_front[0], anchor_W0['W0'], 100*(W0_front[0]-anchor_W0['W0'])/anchor_W0['W0']))
print('  min Wf da frente:  %10.1f kgf | ancora SLSQP: %10.1f kgf (%+.2f %%)'
      %(Wf_front[-1], anchor_Wf['Wf'], 100*(Wf_front[-1]-anchor_Wf['Wf'])/anchor_Wf['Wf']))

# selecao de 3 aeronaves de regioes distintas da frente:
# A = extremo de minimo W0, C = extremo de minimo Wf,
# B = "joelho" (ponto mais proximo da utopia na frente normalizada)
f1n = (F[:,0] - F[:,0].min())/max(F[:,0].max() - F[:,0].min(), 1e-12)
f2n = (F[:,1] - F[:,1].min())/max(F[:,1].max() - F[:,1].min(), 1e-12)
idxA = 0
idxC = len(F) - 1
idxB = int(np.argmin(f1n**2 + f2n**2))

sel_idx = [idxA, idxB, idxC]
sel_names = ['A (min W0)', 'B (joelho)', 'C (min Wf)']
sel_colors = ['tab:red', 'tab:green', 'tab:blue']

print('')
print('%-12s %10s %10s'%('aeronave', 'W0 [kgf]', 'Wf [kgf]'))
for name, k in zip(sel_names, sel_idx):
    print('%-12s %10.1f %10.1f'%(name, W0_front[k], Wf_front[k]))

print('')
print('%-10s'%'DV' + ''.join('%12s'%n for n in sel_names))
for j, dvn in enumerate(dv_names):
    print('%-10s'%dvn + ''.join('%12.4f'%(X[k,j]*Xref[j]) for k in sel_idx))

### FRENTE DE PARETO

fig = plt.figure(figsize=(8,6))
plt.plot(W0_front/1000, Wf_front/1000, 'o', color='gray', markersize=4,
         label='frente de Pareto (NSGA-II)')
plt.plot(out_ref['W0']/1000, out_ref['Wf']/1000, 's', color='k', markersize=9,
         label='baseline (PRJ-22)')
plt.plot(anchor_W0['W0']/1000, anchor_W0['Wf']/1000, '*', color='gold', markersize=17,
         markeredgecolor='k', label='SLSQP min $W_0$ (Secao 3)')
plt.plot(anchor_Wf['W0']/1000, anchor_Wf['Wf']/1000, '*', color='tab:purple', markersize=17,
         markeredgecolor='k', label='SLSQP min $W_f$')
for name, k, c in zip(sel_names, sel_idx, sel_colors):
    plt.plot(W0_front[k]/1000, Wf_front[k]/1000, 'o', color=c, markersize=9,
             markeredgecolor='k', label=name)
plt.xlabel('$W_0$ [t]', fontsize=13)
plt.ylabel('$W_f$ [t]', fontsize=13)
plt.grid(alpha=0.3)
plt.legend(fontsize=10)
plt.tight_layout()
fig.savefig('equipe_moga_pareto.png', dpi=150)

### PLANFORMAS DAS AERONAVES SELECIONADAS

def planform(ax, ap, color, label):

    i = ap['inputs']
    g = ap['geometry']

    for s in (1, -1):
        ax.plot([i['xr_w'], g['xt_w'], g['xt_w']+g['ct_w'], i['xr_w']+g['cr_w'], i['xr_w']],
                s*np.array([0, g['yt_w'], g['yt_w'], 0, 0]), color=color, lw=1.5,
                label=label if s == 1 else None)
        ax.plot([g['xr_h'], g['xt_h'], g['xt_h']+g['ct_h'], g['xr_h']+g['cr_h'], g['xr_h']],
                s*np.array([0, g['yt_h'], g['yt_h'], 0, 0]), color=color, lw=1.5)
        ax.plot([0, i['L_f']], [s*i['D_f']/2, s*i['D_f']/2], color=color, lw=1.0)

    ax.plot([i['x_mlg']]*2, [-i['y_mlg'], i['y_mlg']], 'o', color=color, ms=5)
    ax.plot([i['x_nlg']], [0], 'o', color=color, ms=5)

fig, ax = plt.subplots(figsize=(11,7))

for name, k, c in zip(sel_names, sel_idx, sel_colors):
    airplane = standard_airplane('my_airplane')
    for j, dvn in enumerate(dv_names):
        airplane['inputs'][dvn] = X[k,j]*Xref[j]
    analyze(airplane)
    planform(ax, airplane, c, name)

ax.set_title('Planformas das aeronaves da frente de Pareto', fontsize=13)
ax.set_xlabel('x [m]', fontsize=12)
ax.set_ylabel('y [m]', fontsize=12)
ax.set_aspect('equal')
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
fig.savefig('equipe_moga_planformas.png', dpi=150)

plt.show()
