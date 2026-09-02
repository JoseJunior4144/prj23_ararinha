# Lab 02 — Otimização: resumo de resultados para o relatório

Gerado a partir das execuções dos scripts do repositório (branch `durin`,
com a correção da restrição `mlg_fit` e a nova restrição `ground_clearance`).

As figuras estão organizadas por tópico: `1_monoobj_fokker100/`,
`2_monoobj_equipe/` e `3_multiobj/`.

## Tópico 1 — Mono-objetivo: aeronave default (`lab2_opt_fokker100.py`)

Problema: min MTOW, DVs `AR_w` ∈ [7, 12] e `S_w` ∈ [80, 120] m², restrição
`b_w ≤ 30 m`, SLSQP com diferenças finitas, partida em (7,5; 90 m²).

| | AR_w | S_w [m²] | MTOW [kgf] | b_w [m] |
|---|---|---|---|---|
| Ponto de partida | 7,500 | 90,00 | 43.066,6 | 25,98 |
| Ponto ótimo | 10,522 | 84,88 | 41.879,7 | 29,89 |

- **Melhoria relativa do objetivo:** 2,756 %
- **Chamadas de função:** 33 (`nfev`), 68 execuções do `analyze`, 11 iterações
- **O ótimo é restringido?** Não — é um ótimo **interior**: `b_w = 29,89 m`
  fica aquém do limite de 30 m (g = −3,8·10⁻³, inativa), nenhum bound ativo
  e ‖∇f‖ ≈ 7·10⁻⁶ no ótimo. A condição de otimalidade é ∇f = 0, sem
  contribuição de multiplicadores.

## Tópico 2 — Mono-objetivo: aeronave da equipe (`lab2_opt_equipe_geom.py`)

### Storytelling: da baseline ao ótimo (`equipe_etapas_w0.png`, `equipe_doe_sensibilidade.png`, `equipe_v1v2_delta.png`)

- **v1 → v2** (`equipe_v1v2_delta.png`): o resumo executivo da otimização —
  asa mais esbelta (AR +22,5%), menor (S_w −9,3%), mais enflechada (+5,3%) e
  de raiz mais espessa (+16,3%), com trem mais largo (+26,4%) e mais longo
  (+4,8%), compra **−5,4% de MTOW, −3,3% de peso vazio e −9,5% de
  combustível**, pagando +5,4% de envergadura (58,9 m, ainda código E).
  A v2 está registrada em `Histórico de Versões/standard_airplane v2.py`.

- **Etapas** (`lab2_doe_equipe.py`): baseline PRJ-22 (306,7 t — já inviável:
  trem principal atrás do bordo de fuga, `mlg_fit` = 1,013) → otimização com
  6 DVs (288,3 t — ótimo *aparente*, mas o otimizador empurrou o trem ainda
  mais para trás, `mlg_fit` = 1,106) → formulação final com 8 DVs e 16
  restrições (290,1 t, viável). **As restrições de realismo custam +1,9 t** —
  é o preço de um avião que para em pé.
- **Cortes 1-a-1 no ótimo** (DOE, seguindo o passo-a-passo das aulas):
  `AR_w`, `S_w`, `sweep_w` e `tcr_w` movem o objetivo; `xr_w`, `x_mlg`,
  `z_lg` e `y_mlg` têm curva de W0 **plana** — não são "otimizadas", são
  *posicionadas pelas restrições*. Com 6 restrições ativas, qualquer passo
  numa DV isolada sai do conjunto viável, e a figura mostra **qual restrição
  barra cada direção**: `SM_aft` barra AR/enflechamento baixos, `SM_fwd`
  barra os altos, `tank_excess` barra raiz fina, `alpha_tipback` barra asa
  grande e trem adiantado, `ground_clearance` barra trem curto,
  `mlg_track` barra trem largo.

### 1. Definição do problema e escolha das DVs

min W0 com 8 DVs: `AR_w`, `xr_w`, `S_w`, `sweep_w`, `tcr_w` (asa: arrasto
induzido, de onda e peso estrutural; posição da asa controla SM e tanque) e
`x_mlg`, `z_lg`, `y_mlg` (trem de pouso: viabilizam tipback/tailstrike/
overturn/frações do triquilho ao mesmo tempo em que a asa muda). DVs e
restrições normalizadas pelos valores de referência, como recomendado no
roteiro (essencial para o condicionamento do SLSQP — sem isso a jacobiana
mistura escalas de 0,05 a 390).

### 2. Restrições/objetivos adicionados além do mínimo do roteiro

- `T0req ≤ 1` — empuxo requerido em todas as condições (FAR 25) não pode
  exceder o Tmax dos motores escolhidos;
- `vt_fit ≤ 1` — a raiz da EV deve terminar dentro da fuselagem;
- `mlg_fit ≤ 1` — o trem principal deve ficar sob a asa (à frente do bordo
  de fuga na estação lateral do trem). O bordo de fuga é calculado
  interpolando o LE entre `xr_w` e o `xt_w` do designTool (que já embute o
  termo de 1/4 de corda, pois `sweep_w` é medido a 1/4 de corda);
- `ground_clearance ≥ 0,5 m` — fundo da nacele não pode raspar no solo
  (necessária porque `z_lg` é DV e seu bound superior permitiria nacele
  abaixo do solo);
- Limites do Anexo 14 na interpretação estrita ("up to but not including"):
  `b_w ≤ 64,9 m` e `mlg_track ≤ 13,9 m` (código E); `h_tail ≤ 20 m`
  (FAA ADG V).

### 3. Comparação inicial × otimizado

| DV | inicial | otimizado | variação |
|---|---|---|---|
| AR_w | 8,000 | 9,800 | +22,5 % |
| xr_w [m] | 17,000 | 16,011 | −5,8 % |
| S_w [m²] | 390,00 | 353,67 | −9,3 % |
| sweep_w [rad] | 0,580 | 0,611 | +5,3 % |
| x_mlg [m] | 31,000 | 29,263 | −5,6 % |
| tcr_w | 0,180 | 0,209 | +16,3 % |
| z_lg [m] | −5,750 | −6,025 | +4,8 % |
| y_mlg [m] | 5,500 | 6,950 | +26,4 % |

| Objetivo | inicial | otimizado | variação |
|---|---|---|---|
| W0 [kgf] | 306.745,2 | 290.117,1 | **−5,42 %** |
| Wf [kgf] | 122.743,5 | 111.036,5 | −9,54 % |

Restrições no ótimo (16 no total): 6 ativas — `SM_fwd` (0,30),
`alpha_tipback` (15°), `alpha_tailstrike` (10°), `tank_excess` (0),
`mlg_track` (13,9 m) e `mlg_fit` (1,0). Nenhum bound ativo.
`ground_clearance` = 1,02 m (folgada). Obs.: a baseline do PRJ-22 era
**inviável** em `mlg_fit` (trem em x = 31 m, atrás do bordo de fuga real
≈ 30,6 m) — a otimização corrigiu isso.

### 4. Otimizador

SLSQP (scipy) com gradientes por diferenças finitas, conforme o roteiro.
Adequado porque: problema suave, restrições de desigualdade não lineares,
8 DVs (gradiente escala melhor que livre-de-gradiente), análise barata e
determinística (~7 ms), e o SLSQP trata bounds + restrições nativamente.

### 5. Melhoria relativa

−5,42 % em W0 (−16,6 t) e −9,54 % em Wf (−11,7 t) em relação à aeronave
do PRJ-22.

### 6. Custo

91 chamadas de `objfun` (`nfev`), 182 execuções do `analyze`
(objetivo + restrições em varreduras separadas de diferenças finitas),
10 iterações do SLSQP, ~0,6 s no total.

### 7–9. Figuras

- `equipe_geom_historico.png` — histórico de DVs, objetivo e restrições;
- `equipe_geom_planformas.png` — planta e vista lateral, baseline ×
  otimizado;
- `equipe_geom_3dview.png` — vista 3D da aeronave otimizada.

Física do resultado: AR↑ e S_w↓ reduzem arrasto induzido e molhado; a raiz
mais espessa (`tcr_w` 0,21) barateia a estrutura da asa mais esbelta com
pouca penalidade de onda no modelo (Korn usa t/c médio 0,25·tcr+0,75·tct e
o enflechamento ↑5% compensa); o combustível fecha justo (`tank_excess` = 0)
e o CG dianteiro encosta em SM = 0,30. O trem alarga até o limite de track
de 13,9 m porque a estação mais externa da asa enflechada tem bordo de fuga
mais recuado, aliviando `mlg_fit` e permitindo o trem mais atrás
(tipback/tailstrike ativos amarram `x_mlg` e `z_lg`).

## Tópico 3 — Multiobjetivo W0 × Wf (`lab2_opt_equipe_moga.py`)

min {W0, Wf} com as mesmas 8 DVs e 16 restrições, NSGA-II (pymoo),
**G ≤ 0 no pymoo** (sinal oposto ao scipy). Avaliações que divergem o
ponto-fixo do W0 (população aleatória gera geometrias extremas) são
devolvidas como indivíduos fortemente inviáveis.

### 1–2. Frente de Pareto e parâmetros do algoritmo

- **Indivíduos por geração (pop_size): 80** — **gerações (n_gen): 200** —
  16.000 execuções do `analyze`, ~63 s.
- População inicial: LHS semeado manualmente (o LHS do pymoo não usa o RNG
  global do numpy) + baseline + 2 sementes estritamente viáveis (SLSQP com
  restrições apertadas em 0,5 %). Script reprodutível: duas execuções
  independentes produzem frentes bit a bit idênticas.
- **80 pontos** na frente final (figura `equipe_moga_pareto.png`; dados em
  `equipe_moga_frente.csv`).

### 4. Três aeronaves de regiões distintas da frente

| aeronave | W0 [kgf] | Wf [kgf] | S_w [m²] | tcr_w | AR_w |
|---|---|---|---|---|---|
| A (min W0) | 290.758,2 | 110.994,1 | 358,1 | 0,205 | 9,80 |
| B (joelho) | 291.291,1 | 109.891,3 | 368,8 | 0,196 | 9,81 |
| C (min Wf) | 292.468,2 | 109.385,4 | 376,1 | 0,186 | 9,79 |

De A para C troca-se ~1,7 t de MTOW por ~1,6 t de combustível: asa maior e
mais fina voa mais eficiente (menos arrasto induzido/de onda por área),
mas pesa mais. Em planta as três são quase idênticas
(`equipe_moga_planformas.png`) — o trade-off está em área e espessura,
não na forma.

### Verificação de convergência com a Seção 3 (pergunta 3)

Os extremos da frente devem coincidir com os ótimos mono-objetivo do
mesmo problema: min W0 (SLSQP) = 290.117,1 kgf e min Wf (SLSQP) =
109.117,9 kgf (W0 = 291.961,0 kgf). Uma rodada **ingênua e reprodutível**
(`lab2_moga_convergencia.py`: população inicial aleatória, sem sementes,
pop 60 × 120 gerações) estagna com apenas 2 pontos, **+2,1 % acima da
âncora** e inteiramente dominada pelos ótimos do SLSQP — evidência objetiva
de não-convergência, visível em `equipe_moga_convergencia.png`. A
rodada final (pop 80, 200 gerações, semeada) fecha os extremos em
**+0,22 %** (min W0 = 290.758,2 kgf) e **+0,25 %** (min Wf =
109.385,4 kgf) das âncoras — o resíduo é esperado do NSGA-II em ótimos
"de canto" com várias restrições ativas, que o gradiente resolve com
precisão de máquina.

A frente é estreita porque W0 e Wf são apenas parcialmente conflitantes:
entre os extremos trocam-se ~1,7 t de MTOW por ~1,6 t de combustível —
a aeronave de min Wf tem asa maior (S_w ≈ 376 m²) e raiz mais fina
(tcr ≈ 0,186), mais eficiente em cruzeiro porém estruturalmente mais
pesada.
