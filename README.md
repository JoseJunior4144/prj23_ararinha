# Lab 02 - Otimização (PRJ-23) — Equipe Ararinha

Instituto Tecnológico de Aeronáutica — Programa de Especialização em Engenharia Aeronáutica.
Laboratório de otimização aplicada ao projeto conceitual de aeronaves.

## Como rodar

```
pip install -r requirements.txt
python lab2_opt_fokker100.py       # Topico 1: mono-obj default   (~1,5 min, inclui contorno do espaco de projeto)
python lab2_opt_equipe_geom.py     # Topico 2: mono-obj equipe    (~5 s)
python lab2_doe_equipe.py          # Topico 2: DOE + etapas       (~10 s)
python lab2_opt_equipe_moga.py     # Topico 3: multiobj NSGA-II   (~1-2 min, 16.000 analises)
python lab2_moga_convergencia.py   # Topico 3: evidencia de convergencia (~40 s)
```

Cada script imprime os resultados no terminal e salva as figuras `*.png`
em `Resultados/`. Os scripts são determinísticos (semente fixa no NSGA-II),
então re-rodar reproduz os números do relatório.

## Estrutura

### Exercícios do roteiro
- `lab2_opt_fokker100.py` — **Seção 2**: aeronave default (`fokker100`), min MTOW com
  `AR_w` e `S_w`, restrição `b_w ≤ 30 m`, SLSQP. Resultado: −2,76 % (ótimo interior).
- `lab2_opt_equipe_geom.py` — **Seção 3**: aeronave da equipe (PRJ-22), min W0 com
  8 DVs (asa + trem de pouso) e 16 restrições normalizadas, SLSQP.
  Resultado: **W0 = 290.117 kgf (−5,42 %)**, 6 restrições ativas.
- `lab2_opt_equipe_moga.py` — **Seção 4**: min {W0, Wf} com NSGA-II (pymoo),
  pop. 80 × 200 gerações, população semeada com os ótimos SLSQP (âncoras de
  convergência). Salva a frente em `Resultados/equipe_moga_frente.csv`.
- `lab2_doe_equipe.py` — storytelling do Tópico 2: cortes 1-a-1 no ótimo
  (qual restrição barra cada direção) e comparação de W0 pelas etapas.
- `lab2_moga_convergencia.py` — Tópico 3: rodada ingênua do NSGA-II
  comparada à frente final e às âncoras SLSQP (evidência de convergência).
- `lab2_opt_equipe.py` — versão intermediária da Seção 3 (6 DVs, sem as
  restrições de encaixe do trem), mantida como registro da progressão.

### Ferramenta de projeto
- `designTool/` — módulo de análise conceitual (geometria, aerodinâmica, pesos,
  propulsão, desempenho, balanceamento, trem de pouso). Chamar via `analyze()`.

### Exemplos de aula
- `exemplos_aula/` — `01_uncon_opt_class.py`, `02_con_opt_class.py`,
  `03_multiobj_class.py`, `plot_rosen.py`, `auxmod.py`
- `analyze.py` — demo do designTool (fica na raiz porque importa o pacote)

### Resultados
- `Resultados/` — figuras e dados organizados pelos 3 tópicos do roteiro
  (`1_monoobj_fokker100/`, `2_monoobj_equipe/`, `3_multiobj/`) +
  `lab02_respostas_resumo.md` com as respostas. O relatório final é mantido
  no Overleaf da equipe.

## Notas de modelagem

Duas correções documentadas no apêndice do relatório:
1. O bordo de fuga usado na restrição `mlg_fit` interpola o bordo de ataque
   entre a raiz e o `xt_w` do designTool, pois `sweep_w` é medido a 1/4 de
   corda (`geometry.py`) — a forma anterior subestimava o TE em ~0,45 m e
   custava ~2,9 t de MTOW com a restrição ativa.
2. Restrição `ground_clearance ≥ 0,5 m` adicionada: com `z_lg` como DV, o
   bound superior permitiria a nacele abaixo do solo.

## Dependências

Ver `requirements.txt` (numpy, scipy, matplotlib, pymoo). Testado com
Python 3.12, scipy 1.16, pymoo 0.6.
