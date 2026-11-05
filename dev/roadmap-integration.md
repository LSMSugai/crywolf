# Roadmap de integração sheepdog + crywolf

O plano conjunto está no repositório do sheepdog:
https://github.com/diogoprov/sheepdog/blob/main/dev/roadmap-integration.md
(PR https://github.com/diogoprov/sheepdog/pull/35). Este arquivo só aponta
para lá, para não haver duas versões.

O que toca o crywolf, em ordem de prioridade:

- **Fase 0, contrato de dados:** `fp_read_sheepdog()` lê o diretório
  escrito por `sheepdog::export_occupancy_data()` (esquema
  `sheepdog-occupancy`, v1) e recusa ou avisa quando a validação vem de
  `sampling_mode = "top"`. **Feito em 24/10/2026** (lado do sheepdog no
  PR diogoprov/sheepdog#36): descarta os clipes "top" com aviso
  (`keep_top = TRUE` os mantém), devolve as tabelas no formato de
  `fp_aggregate(by = c("site", "occasion_id"))` e `fp_data()`; teste de
  contrato com a exportação real do sheepdog, que roda no CI quando o
  r-universe publicar o sheepdog 3.2.1.9000.
- **Fase 1, propagação Monte Carlo:** `fp_occu_mc()` ajusta o modelo nos
  históricos de `sheepdog::history_draws()` e junta as posteriores.
  **Feito em 24/10/2026 (crywolf 1.9.0.9000):** `fp_occu_mc()` (mistura
  igual; `summary()` com variância dentro e entre ajustes; com M = 1, draws
  idênticos ao `fp_occu()`) e `fp_draws_data()`, que monta as versões a
  partir de `calibration_draws.parquet` (regra `"sample"` por padrão,
  decisão do Diogo; `"threshold"` como alternativa). Falta: escolher M
  (estabilidade com M = 10, 20, 50) e a comparação na curiango.
  **25/10/2026:** uma primeira rodada (M = 50, `fp = NULL` e `fp ~ 1`) foi
  descartada: o export tinha um artefato (score do minuto repetido em
  cada linha, ~7 detecções por minuto) e a calibração era por janela, não
  por minuto (corrigido em sheepdog#39). Com o export corrigido, a regra
  `"sample"` dá 7,9% dos minutos com detecção (Custom) e `"threshold"`
  1,3%. Os 21,4% das gravações anotadas não servem de referência para o
  total (meses concentrados).
  **Rodada no export corrigido (25/10/2026, `dev/mc-curiango/results.rds`):**
  M = 50, 7,7% dos minutos com detecção por versão. 99 de 100 ajustes
  convergiram (o 50º de `fp ~ temp + rh`, R-hat 1,53, fica fora: valores
  de M = 30). M = 10 já basta. Ocupação, `fp ~ 1`: temperatura −0,07
  (−0,24; 0,11), umidade −0,64, lua 0,15; `fp ~ temp + rh`: temperatura
  −0,15 (−0,30; 0,01), umidade −0,54, lua 0,14, e `fp_temp` = +0,17 (0,09;
  0,24). Limiar 0,4 com `fp ~ 1`: −0,25, −0,78, 0,14. Leitura: 30–70% da
  variância dos efeitos de ocupação vem da calibração (o limiar fixo a
  esconde); o FP volta a crescer com a temperatura, como no limiar 0,2, e
  modelá-lo aproxima o efeito da temperatura do valor do limiar 0,4.
  Ideia: `fp_occu_mc()` poderia avisar ou excluir ajustes que não
  convergem (**feito no 2.0.0:** `max_rhat = 1.05` deixa de fora, com
  aviso, os ajustes cujo maior R-hat passa do limite). **Não** pôr covariáveis da noite em P(TP) como atalho
  (argumento da sessão do sheepdog, 25/10/2026): uma calibração
  discriminativa P(TP | score, temp) absorve a prevalência do canto, que
  também depende da temperatura, e contaria a prevalência duas vezes. O
  caminho com princípio é a Fase 2 (covariáveis nas distribuições dos
  scores f1 e f0, prevalência vinda de ψ e θ); até lá, `fp ~ temp + rh`
  é a correção do lado do crywolf.
  **Ressalva (26/10/2026, PPC das calibrações, sheepdog#41):** a calibração
  do Custom ajusta bem no total, por score e por local, mas erra por mês
  (dez/23 3,5% de TP observados contra 15% previstos; jan/24 0,9% contra
  12%; jun/24 85% contra 74%). Os históricos `"sample"` ganham TP demais nos
  meses de pouca atividade (quentes e úmidos), o que pode explicar a
  atenuação da temperatura e o `fp_temp` positivo sem confusão estruturada
  do classificador. A leitura acima fica como hipótese. Teste possível
  (decisão do Diogo): calibração com estação/mês como análise de
  sensibilidade e nova rodada do `fp_occu_mc()`; se os efeitos voltarem
  para perto dos do limiar 0,4, a atenuação era artefato da calibração.
  **Rodada com calibração por mês (26/10/2026, AWS us-east-2, 27 min,
  `dev/mc-curiango/results-season.rds`; export `curiango-export-season`
  da sessão do sheepdog):** M = 50, 8,8% dos minutos com detecção por
  versão (7,2–10,8%). `fp ~ 1` convergiu nos 50; `fp ~ temp + rh` deixou
  2 de fora (versões 5 e 11, R-hat 1,53). Ocupação, `fp ~ 1`: temperatura
  −0,11 (−0,32; 0,08), umidade −0,55 (−0,82; −0,34), lua 0,02; `fp ~
  temp + rh`: −0,09, −0,49, 0,03, e **`fp_temp` = −0,04 (−0,21; 0,12)**,
  `fp_rh` = −0,14 (−0,32; 0,05). Leitura: (1) o `fp_temp` positivo some
  quando a calibração acompanha o mês, então era erro mensal da calibração,
  não confusão estruturada do classificador; (2) os efeitos de ocupação
  **não** voltam para os do limiar 0,4 (−0,25, −0,78), então a atenuação
  não se explica pelo erro mensal: a diferença para o limiar fixo segue
  sem árbitro (o limiar muda o que conta como detecção); (3) 70–98% da
  variância vem de uma versão para outra (antes 30–70%): por mês há menos
  minutos de teste e a calibração fica mais incerta. Ressalvas: o mês na
  calibração mistura atividade sazonal (prevalência) com o comportamento do
  classificador e reusa os minutos de teste por mês, então é análise de
  sensibilidade, não o modelo principal; o caminho com princípio continua
  sendo a Fase 2.
- **Esquemas 1.1 e 1.2 do export (crywolf 2.1.0, 27/10/2026; lado do
  sheepdog em sheepdog#43 e #44):**
  - 1.1: `fp_read_sheepdog()` aceita qualquer 1.x (2.x é recusado); as
    covariáveis de ocasião viram colunas de `occasions`, e as chaves da
    calibração são lidas como texto. `fp_draws_data()` casa cada chave
    de `schema$calibration$keys` (mês, estação, gravador) pela coluna de
    mesmo nome nas ocasiões; ocasião sem valor de uma chave fica como não
    amostrada. Exports 1.0 calibrados por mês ainda tiram o mês da data.
  - 1.2: `fp_draws_data()` lê os históricos já sorteados pelo sheepdog
    (`occasion_draws.parquet`) quando espécie, regra, semente ou `p_min`
    e ocasiões batem; senão sorteia como antes. Na fixture, o lido é
    idêntico ao calculado.
  - Fixture `tests/testthat/fixtures/sheepdog-export`: cópia do
    `export-1.2` do sheepdog. `dev/check-sheepdog-fixture.R` confere que
    continua idêntica (SHAs de blob pela API do GitHub; `--update` copia).
    Não roda no CI: rodar quando o sheepdog mudar a fixture.
  - Do lado do sheepdog: teste de contrato com o crywolf (item 7,
    sheepdog#46) e fluxo de ponta a ponta do Curiango na vinheta
    `crywolf-integration` com `fp_occu_mc()` e M = 20 (item 6,
    sheepdog#45). O CI do sheepdog não instala o crywolf, então a parte
    do teste que usa o crywolf só roda localmente.
  - Pendente do lado do sheepdog: rodar de novo o modo sazonal do
    `dev/benchmarks/curiango-export.R`, que agora usa
    `calibration_keys = "month"`.
- **Calibração com covariável ambiental ou índice acústico (esquema 1.3,
  aprovado em 28/10/2026):** proposta em
  [`proposal-calibration-covariates.md`](proposal-calibration-covariates.md).
  Com covariável contínua, o sheepdog sorteia as versões no export
  (`occasion_draws.parquet` obrigatório, sem `calibration_draws.parquet`) e
  o crywolf só as lê; análise de sensibilidade 2 × 2 (calibração com e sem
  x, `fp ~ 1` e `fp ~ x`) até um cenário do fpcalib dizer qual é a
  principal; sem pesos 1/`inclusion_prob`. Lado do crywolf depois do artigo
  da MEE.
- **Fase 2, modelo conjunto score–ocupação:** `fp_occu_score()`, função
  nova, com as distribuições dos scores dos verdadeiros e dos falsos
  positivos, o termo de Markov e o desenho de verificação. **Decidido em
  27/10/2026 (Diogo, com a posição da sessão do sheepdog):** entra
  **depois** do `fp_occu_dyn()`, como função própria, sem misturar com o
  dinâmico. Nada do lado do sheepdog pede a Fase 2 antes: os itens da
  integração estão feitos e só a Fase 3 depende dela.
  - Primeiro teste da simulação: a hipótese de que a calibração
    discriminativa P(TP | score) mistura prevalência e comportamento do
    classificador ainda é inferência, não foi testada. Conferir antes de
    usá-la como argumento.
  - Dados: o esquema 1.x atual basta, sem 1.3 nem 2.0. `scores.csv` tem o
    score máximo por ocasião × espécie (censurado em `min_confidence`; par
    ausente = censurado), que é o s_ij do modelo; `detections.csv` tem
    todas as predições, se f1/f0 for por janela; `validation.csv` tem o
    `confidence` de cada clipe (um clipe negativo informa f0 naquele
    score); `recordings_verified.csv` traz as gravações aleatórias com o w
    observado; as covariáveis U_ij entram como `occasion_covariates`. O
    que faltar entra como acréscimo 1.x (por exemplo o score do clipe
    revisado quando difere do da ocasião).
  - Divisão: o sheepdog escreve dados e desenho; f1/f0 e o ajuste ficam
    só no crywolf.
  - **Antes de começar o `fp_occu_score()`, ler no sheepdog
    `dev/dag/NOTAS-fase2.md`** (sessão de 02–04/11/2026). Ele resume a
    leitura completa de Rhinehart et al. (2022) e Ogawa et al. (2025):
    de onde vem cada parte do modelo, as lições de identificação e as
    decisões em aberto. No mesmo diretório estão a DAG de referência, no
    estilo da DAG deste pacote (`fase2_dag.dot`/`.png`), e a tabela
    comparativa (`comparacao_modelos.md`). O modelo, já na notação da
    vinheta do crywolf (i, k, **x**, **d**, **g**, gap Δ_s), está no
    roadmap do sheepdog, na seção "O modelo (proposta, a validar)".
  - **Correção (04/11/2026):** o score máximo de `scores.csv` não pode
    entrar como se seguisse f1/f0 da janela. O máximo de n janelas tem
    distribuição própria, que depende de `n_windows`; ver a seção
    "Distribuição do máximo" no roadmap do sheepdog.
  - A hipótese de que a calibração mistura prevalência e classificador já
    tem evidência empírica no Curiango (sheepdog#67), mas a simulação da
    Fase 2 ainda deve conferir.
- **Fase 3:** `fp_review_priority()` para a camada dirigida da validação.
- **Fases 4 e 5:** atividade vocal e comunidade, alinhadas com a Fase A do
  multiespécies do [roadmap v2](roadmap-v2.md). A confusão entre espécies
  continua fora do plano.

Nada disso muda a API congelada no 2.0.0: tudo entra em funções novas.
