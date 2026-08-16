# Bakta Distributed Annotation

Distributed bacterial genome annotation with Bakta and auditable protein prioritization for functional refinement.

## Visão geral

Este repositório contém dois componentes externos ao Bakta usados em um fluxo didático e reprodutível de anotação genômica bacteriana:

1. **`anot00_bakta_01_fleet.py`** — organiza amostras, cria uma fila compartilhada coordinator–worker e executa o Bakta de forma distribuída.
2. **`anot00_bakta_02_selecionados.py`** — lê as saídas do Bakta, verifica a consistência entre arquivos e prioriza proteínas para refinamento funcional posterior.

> O Bakta realiza a anotação genômica. A distribuição de tarefas, a auditoria e a classificação em prioridades são funções deste pipeline externo e não fazem parte das funções nativas do Bakta.

## Estrutura

```text
bakta-distributed-annotation/
├── README.md
├── .gitignore
├── scripts/
│   ├── anot00_bakta_01_fleet.py
│   └── anot00_bakta_02_selecionados.py
└── docs/
    └── tutorial.md
```

## Requisitos principais

- Linux
- Python 3
- Bakta
- banco de dados compatível com a versão instalada do Bakta
- Conda ou Mamba é recomendado para instalação do ambiente

O segundo script utiliza somente a biblioteca padrão do Python e requer Python 3.9 ou superior.

## Uso rápido

Consulte primeiro a ajuda dos scripts:

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

Fluxo básico:

```text
Genomas do NCBI
      ↓
anot00_bakta_01_fleet.py init
      ↓
fila de jobs
      ↓
anot00_bakta_01_fleet.py work
      ↓
Bakta
      ↓
output00_bakta
      ↓
anot00_bakta_02_selecionados.py
      ↓
proteínas priorizadas para análises posteriores
```

O tutorial de instalação e execução está em [`docs/tutorial.md`](docs/tutorial.md).

## Entradas esperadas

O primeiro script trata cada subdiretório imediato da raiz de entrada como uma amostra e procura recursivamente:

- `*_genomic.fna`
- `*_genomic.gbff`

FASTAs contendo apenas CDS ou RNA do NCBI são ignorados explicitamente.

## Saídas principais

A etapa de seleção classifica CDS elegíveis em:

- `high`
- `medium`
- `low`

Pseudogenes e sORFs são, por padrão, tratados separadamente. O arquivo `<amostra>.priority_high.faa` foi concebido como entrada da etapa seguinte de busca curada.

## Reprodutibilidade

O fluxo inclui verificações de identidade entre JSON e FAA, hashes SHA-256, escrita atômica e manifestos por amostra para aumentar a rastreabilidade das execuções.

## Bakta

Bakta é desenvolvido por Schwengers e colaboradores. Consulte o projeto oficial e cite o artigo original ao utilizar o software em trabalhos científicos.

## Licença

A licença deste repositório ainda será definida. Até que um arquivo `LICENSE` seja adicionado, não se presume uma licença de reutilização para estes scripts.
