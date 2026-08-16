# Tutorial básico de instalação e execução

Este tutorial apresenta os comandos essenciais para instalar o Bakta e executar os dois scripts deste repositório.

## 1. Instalar o Bakta

Com Conda:

```bash
conda create -n bakta_env -c conda-forge -c bioconda bakta
conda activate bakta_env
```

Com Mamba:

```bash
mamba create -n bakta_env -c conda-forge -c bioconda bakta
mamba activate bakta_env
```

Verifique a instalação:

```bash
bakta --version
bakta --help
```

## 2. Instalar o banco de dados

Consulte as versões disponíveis:

```bash
bakta_db list
```

Crie uma pasta para o banco:

```bash
mkdir -p ~/bakta_db
```

Banco completo:

```bash
bakta_db download --output ~/bakta_db --type full
```

Banco reduzido:

```bash
bakta_db download --output ~/bakta_db --type light
```

O primeiro script resolve o caminho do banco nesta ordem:

```text
--db
↓
variável BAKTA_DB
↓
caminho padrão definido no script
```

Também é possível definir:

```bash
export BAKTA_DB=/caminho/para/bakta/db
```

## 3. Consultar a ajuda dos scripts

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

## 4. Organização das amostras

Por padrão, o primeiro script espera uma raiz contendo uma subpasta por amostra. Dentro de cada amostra ele procura recursivamente:

```text
*_genomic.fna
*_genomic.gbff
```

FASTAs contendo apenas CDS ou RNA são ignorados.

Exemplo:

```text
etapa00_genomas_ncbi/
├── amostra_01/
│   ├── arquivo_genomic.fna
│   └── arquivo_genomic.gbff
├── amostra_02/
└── amostra_03/
```

Os caminhos padrão foram definidos para o ambiente original do pipeline. Em outro computador, use `--input-root`, `--output-root` e `--db` para indicar seus próprios caminhos.

## 5. Criar a fila de tarefas

Usando os caminhos padrão:

```bash
python3 scripts/anot00_bakta_01_fleet.py init
```

Usando caminhos próprios:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /caminho/etapa00_genomas_ncbi \
  --output-root /caminho/output00_bakta \
  init
```

## 6. Verificar a fila

```bash
python3 scripts/anot00_bakta_01_fleet.py status
```

Monitoramento contínuo:

```bash
python3 scripts/anot00_bakta_01_fleet.py status --watch
```

Com caminhos próprios:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /caminho/output00_bakta \
  status --watch
```

## 7. Testar uma única amostra

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --cpus 4 \
  --once
```

Com banco explícito:

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --cpus 4 \
  --db /caminho/para/bakta/db \
  --once
```

## 8. Processar a fila até o fim

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --cpus 4 \
  --skip-plot \
  --until-empty
```

Em outro computador, use outro identificador:

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-02 \
  --cpus 4 \
  --skip-plot \
  --until-empty
```

Todos os workers devem enxergar a mesma fila compartilhada, mas cada worker deve ter acesso local a uma instalação funcional do Bakta e ao banco correspondente.

## 9. Parâmetros importantes da etapa Bakta

Algumas opções oferecidas pelo script são:

```text
--cpus N
--db CAMINHO
--complete
--compliant
--keep-contig-headers
--skip-plot
--gram +|-|?
--translation-table 11|4|25
--tmp-dir CAMINHO
--once
--until-empty
```

A opção `--complete` deve ser usada somente quando todas as sequências presentes naquele FASTA forem replicons completos.

## 10. Executar a seleção posterior

Depois que as anotações do Bakta estiverem concluídas:

```bash
python3 scripts/anot00_bakta_02_selecionados.py --dry-run
```

Se a conferência estiver correta:

```bash
python3 scripts/anot00_bakta_02_selecionados.py
```

Com caminhos próprios:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /caminho/output00_bakta \
  --output-root /caminho/output00_bakta_selecao
```

## 11. Processar apenas uma amostra

Em lote, por nome:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --sample amostra_01
```

Mais de uma amostra:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --sample amostra_01 \
  --sample amostra_02
```

Modo direto:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-dir /caminho/output00_bakta/amostra_01 \
  --outdir /caminho/output00_bakta_selecao/amostra_01
```

## 12. Política padrão de seleção

Os valores padrão definidos no script são:

```text
proteína curta:           < 90 aminoácidos  → flag descritiva
identidade fraca:         < 90%             → prioridade baixa quando aplicável
cobertura da query fraca: < 80%             → prioridade baixa quando aplicável
cobertura do alvo fraca:  < 80%             → prioridade baixa quando aplicável
```

Pseudogenes e sORFs são separados por padrão. Produtos hipotéticos ou não caracterizados recebem prioridade alta; produtos genéricos, de domínio, família, DUF ou UPF recebem prioridade intermediária; descrições incertas ou inferências abaixo dos limiares podem receber prioridade baixa.

## 13. Fluxo mínimo

```text
Genomas do NCBI
      ↓
anot00_bakta_01_fleet.py init
      ↓
fila compartilhada
      ↓
anot00_bakta_01_fleet.py work
      ↓
Bakta
      ↓
output00_bakta
      ↓
anot00_bakta_02_selecionados.py --dry-run
      ↓
anot00_bakta_02_selecionados.py
      ↓
output00_bakta_selecao
      ↓
proteínas priorizadas para refinamento funcional
```

## Observação conceitual

O Bakta realiza a anotação genômica. A distribuição de tarefas, a auditoria dos arquivos e a classificação das proteínas em prioridades são funções dos scripts deste repositório e não fazem parte das funções nativas do Bakta.
