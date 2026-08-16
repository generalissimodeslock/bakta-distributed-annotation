# Tutorial básico de instalação e execução

Este tutorial apresenta os comandos essenciais para instalar o Bakta e executar os dois scripts deste repositório.

## 1. Instalar o ambiente

Com Mamba:

```bash
mamba env create -f environment.yml
mamba activate bakta_env
```

Com Conda:

```bash
conda env create -f environment.yml
conda activate bakta_env
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

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /caminho/etapa00_genomas_ncbi \
  --output-root /caminho/output00_bakta \
  init
```

A fila é criada em `output00_bakta/bakta_queue/`, com os estados `pending`, `running`, `done`, `failed` e a pasta `logs`.

## 6. Verificar a fila

Uma consulta simples:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /caminho/output00_bakta \
  status
```

Monitoramento contínuo:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /caminho/output00_bakta \
  status --watch
```

## 7. Testar uma única amostra

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /caminho/output00_bakta \
  work \
  --id worker-01 \
  --cpus 4 \
  --db /caminho/para/bakta/db \
  --once
```

A opção `--once` executa no máximo um job e encerra, sendo útil para validar um novo worker.

## 8. Processar a fila até o fim

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root /caminho/output00_bakta \
  work \
  --id worker-01 \
  --cpus 4 \
  --db /caminho/para/bakta/db \
  --skip-plot \
  --until-empty
```

Em outro computador, use outro identificador, por exemplo `worker-02`.

Todos os workers devem enxergar a mesma fila compartilhada e os mesmos diretórios de entrada/saída. Cada worker deve ter uma instalação funcional do Bakta e acesso a um banco compatível; quando possível, uma cópia local idêntica do banco é preferível para reduzir leitura pela rede.

## 9. Parâmetros importantes da etapa Bakta

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

A opção `--complete` deve ser usada somente quando **todas as sequências presentes no FASTA forem replicons completos**.

## 10. Notificação opcional por e-mail

As configurações devem ser fornecidas explicitamente:

```bash
export BAKTA_EMAIL_TO="destinatario@example.org"
export BAKTA_SMTP_USER="remetente@example.org"
export BAKTA_GMAIL_APP_PASSWORD="senha-de-aplicativo"
```

Depois execute:

```bash
python3 scripts/anot00_bakta_01_fleet.py work \
  --id worker-01 \
  --notify-email \
  --until-empty
```

Nunca versione senhas, arquivos `.env` ou outros segredos no repositório.

## 11. Conferir a seleção posterior

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /caminho/output00_bakta \
  --output-root /caminho/output00_bakta_selecao \
  --dry-run
```

O modo `--dry-run` verifica a descoberta das amostras e dos arquivos necessários sem gerar novas saídas.

## 12. Executar a seleção posterior

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /caminho/output00_bakta \
  --output-root /caminho/output00_bakta_selecao
```

## 13. Processar apenas uma amostra

Em lote, por nome:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root /caminho/output00_bakta \
  --output-root /caminho/output00_bakta_selecao \
  --sample amostra_01
```

Modo direto:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-dir /caminho/output00_bakta/amostra_01 \
  --outdir /caminho/output00_bakta_selecao/amostra_01
```

## 14. Política padrão de seleção

```text
proteína curta:           < 90 aminoácidos  → flag descritiva
identidade fraca:         < 90%             → prioridade baixa quando aplicável
cobertura da query fraca: < 80%             → prioridade baixa quando aplicável
cobertura do alvo fraca:  < 80%             → prioridade baixa quando aplicável
```

Pseudogenes e sORFs são separados por padrão. Produtos hipotéticos ou não caracterizados recebem prioridade alta; produtos genéricos, de domínio, família, DUF ou UPF recebem prioridade intermediária; descrições incertas ou inferências abaixo dos limiares podem receber prioridade baixa.

Flags como proteína curta ou ausência de EC, COG e GO servem para auditoria e não representam, por si só, erro de anotação.

## 15. Saídas principais da seleção

```text
<amostra>.targets_refined.tsv
<amostra>.selection_audit.tsv
<amostra>.priority_high.faa
<amostra>.priority_medium.faa
<amostra>.priority_low.faa
<amostra>.sorf_separate.faa
<amostra>.pseudogene_separate.faa
<amostra>.target_report.tsv
<amostra>.target_manifest.json
```

## 16. Fluxo mínimo

```text
Genomas do NCBI
      ↓
inicialização da fila
      ↓
execução distribuída do Bakta
      ↓
saídas do Bakta
      ↓
dry-run da seleção
      ↓
auditoria e seleção
      ↓
FASTAs por prioridade + tabelas de auditoria + manifesto
```

## Observação conceitual

O Bakta realiza a anotação genômica. A distribuição de tarefas, a auditoria dos arquivos e a classificação das proteínas em prioridades são funções dos scripts deste repositório e não fazem parte das funções nativas do Bakta.
