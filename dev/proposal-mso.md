# Proposta: modelo multiespécies com falsos positivos (`fp_occu_ms()`, Fase A)

Escrita em 02/11/2026 pela sessão do crywolf, a partir da Fase A do
`roadmap-v2.md` e da Fase 5 do `roadmap-integration.md` do sheepdog
(sheepdog#74). **Aprovada pelo Diogo em 02/11/2026, como escrita**: as
quatro perguntas do fim ficam com a opção proposta (lista fixa; sem LKJ no
começo; gravador comum; `ver` por espécie no PR 3).

## Resumo

1. Função nova `fp_occu_ms()`. A API do `fp_occu()` não muda.
2. Cada espécie tem ocupação, detecção verdadeira e falso positivo
   próprios, com os coeficientes sorteados de distribuições da comunidade
   (hierárquico, parametrização não centrada). O falso positivo de uma
   espécie é independente do das outras: o classificador erra ao marcar a
   espécie, não ao trocá-la por outra.
3. A lista de espécies é fixa (a do export ou a dos dados). Sem aumento de
   dados para espécies nunca detectadas, por enquanto: a riqueza estimada é
   a das espécies da lista.
4. Caminho incremental que reaproveita a verossimilhança já testada do
   `fp_occu()`, como o `fp_occu_dyn()` fez com as temporadas.
5. Cinco PRs. Validação por simulação e SBC na AWS, com a receita de
   `dev/dyn-sim/aws/`.

## 1. Dados

`fp_data(species = "species")`: coluna de espécie em `auto` e em
`verified`, formato longo (uma linha por ocasião e espécie). Cada ocasião
gravada aparece para todas as espécies da lista, com `y = 0` onde a espécie
não foi detectada. Sem isso não há não-detecções.

- `print()` e os avisos de desenho por espécie: detecções, sítios com
  detecção, ocasiões verificadas e verificações positivas por espécie.
- Aviso quando uma espécie não tem verificação nenhuma: o falso positivo
  dela só é identificado pela comunidade e pelo prior.
- `fp_read_sheepdog(species = NULL)` passa a manter todas as espécies sem
  o aviso atual quando o destino é o multiespécies; `fp_draws_data()` lê
  as versões de todas as espécies juntas (hoje para em
  `fp_draws_data.R:120`).
- Teste sobre a fixture `export-1.3-ms` (já copiada,
  `tests/testthat/fixtures/sheepdog-export-1.3-ms`): três espécies, uma
  rara (só em ARA e BAI).

## 2. Verossimilhança

O `fp_occu_dyn()` expande os sítios em unidades sítio × temporada e passa
as unidades para a verossimilhança do `fp_occu()`
(`fp_emissions.stanfunctions`). O multiespécies faz o mesmo com unidades
sítio × espécie, e z é marginalizado por unidade.

A verossimilhança já aceita, por unidade, um deslocamento do intercepto da
detecção (`d_det`, do `random_det`) e do falso positivo (`e_fp`, do
`random_fp`). Os interceptos por espécie entram por aí, sem mexer nas
emissões:

- `d_det[unidade] = alpha0[espécie]` e `e_fp[unidade] = gamma0[espécie]`,
  com o gap por espécie (`alpha0_s - gamma0_s > 0`), como no `fp_occu()`;
- a ocupação fica fora das emissões, então interceptos e declives por
  espécie são diretos.

Assim, o primeiro modelo (PR 2) tem interceptos por espécie em occ, det e
fp, declives por espécie em occ e declives comuns em det e fp. Declives por
espécie em det e fp exigem mudar as emissões (coeficientes por linha) e
ficam para o PR 3, com teste de igualdade contra o PR 2.

Pronto quando: com uma espécie, a log-densidade é igual à do `fp_occu()`
(diferença ~1e-12); com várias e sem hierarquia, igual à soma de ajustes
separados por espécie.

## 3. Hierarquia e priors

- Para cada coeficiente por espécie: média e desvio da comunidade,
  `beta_s = mu + sigma * z_s`, `z_s ~ Normal(0, 1)`.
- Priors novos em `fp_priors()`: média da comunidade com os priors atuais
  de intercepto e declive; desvio `sigma ~ Normal+(0, s)`, com `s` a
  definir por simulação de prior preditiva.
- Interceptos de occ, det e fp sem correlação entre si no começo. A
  correlação (LKJ) entre ocupação e detecção é comum em MSOM e fica como
  opção depois da validação.
- `random_recorder`: efeito do gravador comum a todas as espécies (o
  gravador ruim é ruim para todas).
- `markov`, `random_fp`, `random_det` e o dinâmico ficam fora da Fase A;
  `fp_occu_ms()` dá erro se pedidos.

## 4. Saídas

- `summary()` com as médias e desvios da comunidade e os coeficientes por
  espécie (`occ_(Intercept)[espécie]`, etc.), nomes em `crywolf-output`.
- Quantidades derivadas no Stan: riqueza por sítio (soma dos z amostrados,
  pelo condicional de cada unidade), `p10` e `p11` por espécie.
- `fp_sites()` por sítio e espécie, `fp_ppc()` por sítio e espécie,
  `simulate()` e um simulador `simulate_fp_ms()`.

## 5. Validação

- Recuperação com 5, 20 e 50 espécies; cenários com uma espécie rara e com
  espécies sem verificação.
- SBC (1.000 réplicas) sobre médias e desvios da comunidade.
- O que se quer medir: quanto a comunidade empresta à espécie rara, no
  viés e na cobertura do p10 e da ocupação dela, contra o modelo de uma
  espécie só.
- Custo e mixing: tempo e R-hat com 20 espécies num conjunto do tamanho
  da curiango; documentar `threads_per_chain`. Estimativa a medir no
  piloto: o tempo cresce mais ou menos linearmente com o número de
  espécies.

## 6. Divisão em PRs

| PR | Conteúdo | Pronto quando |
|---|---|---|
| 1 | Dados: `fp_data(species = )`, `print()` e avisos por espécie, leitura do export multiespécies, fixture | testes sobre `export-1.3-ms`; avisos por espécie |
| 2 | Stan e `fp_occu_ms()` com interceptos por espécie (via unidades) e hierarquia | igualdade com `fp_occu()` e com ajustes separados; ajuste curto recupera os parâmetros |
| 3 | Declives por espécie em det e fp (emissões por linha) | igualdade com o PR 2 quando os declives são comuns |
| 4 | Saídas: riqueza, `fp_sites()`, `fp_ppc()`, `simulate()`, `simulate_fp_ms()` | testes como os do dinâmico |
| 5 | Simulações e SBC na AWS; vinheta | tabelas em `dev/ms-sim/`; vinheta com dados simulados e, depois do depósito, o FrogNet |

## 7. Dados reais

Segundo o `roadmap-integration.md` do sheepdog (decidido em 02/11/2026), o
caso real é o FrogNet (DiverSons, 7 sítios, 2021–2025), num export fora dos
repositórios dos pacotes; o VicFrogNET foi descartado por não publicar a
lista de gravações. A vinheta usa dados simulados até o depósito do
FrogNet.

## Perguntas (decididas em 02/11/2026 pela opção proposta)

1. Lista fixa de espécies, sem aumento de dados para espécies nunca
   detectadas. Com um classificador, a lista é a das classes do modelo, e
   uma espécie nunca detectada ainda entra com y = 0 em toda parte. Serve?
2. Interceptos sem correlação entre submodelos no começo, com a correlação
   como opção depois. Ou já com LKJ?
3. Efeito do gravador comum a todas as espécies. Ou por espécie?
4. Verificação (`ver`) com intercepto por espécie hierárquico, como os
   outros? A alternativa é um r11 comum, mais simples, mas o revisor pode
   errar mais em umas espécies que em outras. A verossimilhança não tem
   hoje deslocamento por unidade na verificação, então o intercepto por
   espécie em `ver` entraria no PR 3, com os declives por espécie.
