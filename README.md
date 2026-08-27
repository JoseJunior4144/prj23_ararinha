# Lab 02 - Otimização 1 (PRJ-23)

Instituto Tecnológico de Aeronáutica — Programa de Especialização em Engenharia Aeronáutica.
Laboratório de otimização aplicada ao projeto conceitual de aeronaves.

## Estrutura

### Exemplos de aula
- `01_uncon_opt_class.py` — otimização irrestrita
- `02_con_opt_class.py` — otimização com restrições
- `03_multiobj_class.py` — otimização multi-objetivo (NSGA-II / pymoo)
- `plot_rosen.py` — visualização da função de Rosenbrock

### Exercícios
- `lab2_opt_fokker100.py` — otimização do Fokker 100
- `lab2_opt_equipe.py` — otimização da aeronave da equipe (etapa 1)
- `lab2_opt_equipe_geom.py` — otimização geométrica da aeronave da equipe (etapa 2)
- `analyze.py`, `auxmod.py` — funções auxiliares

### Ferramenta de projeto
- `designTool/` — módulo de análise conceitual de aeronaves (aerodinâmica, peso, propulsão, desempenho, estabilidade)

### Resultados
Figuras `*.png` geradas pelos scripts (históricos de convergência, espaço de projeto, planformas, vista 3D).

## Dependências

```
numpy
scipy
matplotlib
pymoo>=0.6.0
```
