# Tutorial de instalação e uso

Este tutorial apresenta o fluxo completo para instalar e utilizar os dois scripts do repositório **Bakta Distributed Annotation**. O objetivo é permitir que uma pessoa que esteja chegando ao projeto consiga reproduzir o fluxo em um único computador ou em vários workers conectados a uma fila compartilhada.

O repositório contém duas etapas principais:

1. `anot00_bakta_01_fleet.py` — descobre as amostras, cria a fila de tarefas, executa o Bakta e acompanha os workers.
2. `anot00_bakta_02_selecionados.py` — audita as saídas do Bakta e prioriza proteínas para refinamento funcional posterior.

> O Bakta realiza a anotação genômica. A fila distribuída, a auditoria entre arquivos e a classificação das proteínas são funções externas implementadas pelos scripts deste repositório.

---

## 1. Requisitos

Antes de começar, o computador deve possuir:

- Linux;
- Git;
- Conda ou Mamba;
- acesso à internet para instalar o ambiente e baixar o banco do Bakta;
- espaço em disco suficiente para o banco do Bakta e para as saídas das análises.

O ambiente fornecido pelo repositório utiliza:

```text
Python >= 3.9
Bakta 1.12.x
```

O banco do Bakta **não é incluído no repositório**, pois é grande e possui versionamento próprio.

---

## 2. Baixar o repositório

No terminal:

```bash
git clone https://github.com/generalissimodeslock/bakta-distributed-annotation.git
cd bakta-distributed-annotation
```

A partir deste ponto, os comandos deste tutorial pressupõem que o terminal está dentro da pasta do repositório.

A estrutura principal será semelhante a:

```text
bakta-distributed-annotation/
├── README.md
├── LICENSE
├── environment.yml
├── docs/
├── scripts/
│   ├── anot00_bakta_01_fleet.py
│   └── anot00_bakta_02_selecionados.py
└── tests/
    └── test_smoke.py
```

---

## 3. Criar o ambiente do programa

### Com Mamba

```bash
mamba env create -f environment.yml
mamba activate bakta_env
```

### Com Conda

```bash
conda env create -f environment.yml
conda activate bakta_env
```

Para confirmar que o ambiente foi ativado corretamente:

```bash
python --version
bakta --version
```

Também é possível consultar a ajuda do Bakta:

```bash
bakta --help
```

---

## 4. Verificar os scripts

Consulte as versões:

```bash
python3 scripts/anot00_bakta_01_fleet.py --version
python3 scripts/anot00_bakta_02_selecionados.py --version
```

Consulte todas as opções disponíveis:

```bash
python3 scripts/anot00_bakta_01_fleet.py --help
python3 scripts/anot00_bakta_02_selecionados.py --help
```

Esses comandos não executam análises; apenas confirmam que os scripts podem ser carregados pelo Python.

---

## 5. Instalar o banco de dados do Bakta

Consulte as versões de banco compatíveis:

```bash
bakta_db list
```

Crie uma pasta para os bancos:

```bash
mkdir -p ~/bakta_db
```

Para baixar o banco completo:

```bash
bakta_db download --output ~/bakta_db --type full
```

Para utilizar a versão reduzida:

```bash
bakta_db download --output ~/bakta_db --type light
```

Depois do download, identifique o diretório real do banco criado pelo `bakta_db` e utilize esse caminho nos comandos abaixo.

O script `anot00_bakta_01_fleet.py` procura o banco nesta ordem:

```text
1. --db
2. variável de ambiente BAKTA_DB
3. caminho padrão presente no script
```

Para tornar o caminho disponível no terminal atual:

```bash
export BAKTA_DB=/caminho/real/para/o/banco/bakta
```

Para uso em outro computador, é recomendado informar `--db` explicitamente ou definir `BAKTA_DB` nesse worker.

---

## 6. Preparar as amostras

O script considera cada subpasta imediata da raiz de entrada como uma amostra.

Exemplo:

```text
etapa00_genomas_ncbi/
├── Acinetobacter_schindleri_HZE30-1/
│   ├── GCF_xxxxx_genomic.fna
│   └── GCF_xxxxx_genomic.gbff
├── Delftia_tsuruhatensis_CM13/
│   ├── GCF_yyyyy_genomic.fna
│   └── GCF_yyyyy_genomic.gbff
└── outra_amostra/
    ├── arquivo_genomic.fna
    └── arquivo_genomic.gbff
```

Dentro de cada pasta, o script procura recursivamente arquivos do tipo:

```text
*_genomic.fna
*_genomic.gbff
```

FASTAs do NCBI contendo apenas CDS ou RNA são ignorados, como:

```text
*_cds_from_genomic.fna
*_rna_from_genomic.fna
```

O arquivo `GBFF` é utilizado para recuperar metadados, enquanto o `FNA` contém a montagem que será enviada ao Bakta.

---

## 7. Definir os diretórios de trabalho

Os scripts ainda preservam caminhos padrão do ambiente em que o pipeline foi desenvolvido. Para uso em outro computador, é mais seguro indicar seus próprios caminhos.

Neste tutorial, considere:

```text
/caminho/etapa00_genomas_ncbi
/caminho/output00_bakta
/caminho/output00_bakta_selecao
```

Substitua `/caminho/` pelo local real utilizado em sua máquina.

Você também pode criar variáveis de conveniência no shell:

```bash
INPUT_ROOT=/caminho/etapa00_genomas_ncbi
BAKTA_OUTPUT=/caminho/output00_bakta
SELECTION_OUTPUT=/caminho/output00_bakta_selecao
```

Essas três variáveis são apenas atalhos do shell; os scripts não dependem delas.

---

# Parte I — Execução do Bakta

## 8. Criar a fila de tarefas

A fila deve ser criada antes de iniciar os workers.

Usando os caminhos diretamente:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root /caminho/etapa00_genomas_ncbi \
  --output-root /caminho/output00_bakta \
  init
```

Ou, se você criou as variáveis de conveniência:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init
```

A fila será criada dentro de:

```text
output00_bakta/
└── bakta_queue/
    ├── pending/
    ├── running/
    ├── done/
    ├── failed/
    └── logs/
```

Cada amostra válida gera um pequeno arquivo JSON de trabalho, ou `job`.

### Recriar a fila

Existe a opção:

```bash
--overwrite
```

Por exemplo:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init --overwrite
```

Essa opção limpa os estados da fila antes de recriá-la. **Não deve ser usada rotineiramente**; utilize apenas quando houver intenção explícita de reconstruir a fila.

---

## 9. Consultar o estado da fila

Uma consulta única:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status
```

Monitoramento contínuo:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status --watch
```

Por padrão, o modo `--watch` é atualizado a cada cinco segundos.

Outro intervalo pode ser definido com:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status --watch --interval 10
```

O painel informa quantos jobs estão em:

```text
pending
running
done
failed
```

Também mostra os workers atualmente ativos, a amostra em processamento e uma estimativa de progresso obtida das mensagens do Bakta.

---

## 10. Testar um único job

Antes de liberar toda a fila, execute somente uma amostra:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --once
```

A opção:

```text
--once
```

faz o worker reivindicar no máximo um job e encerrar depois da execução.

Esse é o modo recomendado para validar uma nova instalação antes de processar todas as amostras.

---

## 11. Executar a fila inteira em um único computador

O uso distribuído **não é obrigatório**. Um único computador pode processar toda a fila:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

A opção:

```text
--until-empty
```

mantém o worker ativo até que não existam mais jobs em `pending` ou `running`.

---

## 12. Executar em vários computadores

Para utilizar vários workers, todos devem enxergar:

- a mesma raiz de entrada;
- a mesma raiz `output00_bakta`;
- a mesma fila `bakta_queue`.

Cada computador executa o mesmo comando, mudando apenas o identificador do worker e, quando necessário, o caminho local do banco.

### Worker 01

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

### Worker 02

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-02 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --until-empty
```

Os identificadores devem ser diferentes. Eles podem ser nomes como `worker-01`, `worker-02` ou endereços IP.

O banco do Bakta não precisa necessariamente estar no mesmo caminho físico em todos os computadores: cada worker pode apontar `--db` para sua própria cópia local compatível.

Para bancos de leitura intensa, cópias locais idênticas tendem a ser preferíveis a uma única cópia servida pela rede.

---

## 13. Principais opções da etapa Bakta

As opções mais úteis do modo `work` são:

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
--no-dashboard
--dashboard-interval SEGUNDOS
```

### `--complete`

Utilize somente quando **todas as sequências presentes no FASTA forem replicons completos**.

Não utilize essa opção indiscriminadamente em montagens WGS fragmentadas em contigs.

### `--skip-plot`

Desativa a produção dos gráficos PNG/SVG do Bakta. É útil quando o objetivo principal são os arquivos de anotação e quando se deseja reduzir saídas desnecessárias.

### `--gram`

Valores aceitos:

```text
+
-
?
```

O padrão é `?`.

### `--translation-table`

Valores aceitos:

```text
11
4
25
```

O padrão é `11`.

### `--tmp-dir`

Permite direcionar os arquivos temporários para um disco local rápido:

```bash
--tmp-dir /caminho/temporario
```

---

## 14. Painel do worker

Durante a execução, o worker possui um painel de terminal próprio.

Para desativá-lo:

```bash
--no-dashboard
```

Para mudar a frequência de atualização:

```bash
--dashboard-interval 2
```

O valor é informado em segundos.

---

## 15. Notificação opcional por e-mail

O repositório público não contém endereço pessoal nem senha de e-mail embutidos no código.

Para habilitar a notificação, defina as informações explicitamente.

Exemplo com as variáveis padrão:

```bash
export BAKTA_EMAIL_TO="destinatario@example.org"
export BAKTA_SMTP_USER="remetente@example.org"
export BAKTA_GMAIL_APP_PASSWORD="senha-de-aplicativo"
```

Depois, escolha **um** worker para ser responsável pela notificação:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work \
  --id worker-01 \
  --cpus 4 \
  --db "$BAKTA_DB" \
  --skip-plot \
  --notify-email \
  --until-empty
```

Também é possível fornecer o destinatário e o usuário SMTP diretamente:

```text
--email-to
--smtp-user
--smtp-host
--smtp-port
--smtp-password-env
```

A senha deve permanecer em uma variável de ambiente e **nunca deve ser adicionada ao GitHub**.

O `.gitignore` do repositório já ignora arquivos `.env` e diretórios de segredos, mas isso não substitui o cuidado do usuário.

---

## 16. Verificar o término da anotação

Ao terminar, cada job será movido para:

```text
done/
```

ou:

```text
failed/
```

Os logs ficam em:

```text
output00_bakta/bakta_queue/logs/
```

Para uma conferência final:

```bash
python3 scripts/anot00_bakta_01_fleet.py \
  --output-root "$BAKTA_OUTPUT" \
  status
```

O script só considera uma execução bem-sucedida quando o Bakta retorna código de saída adequado e os arquivos obrigatórios esperados estão presentes.

---

# Parte II — Auditoria e seleção das proteínas

## 17. O que a segunda etapa utiliza

O script `anot00_bakta_02_selecionados.py` utiliza principalmente:

```text
<amostra>.json
<amostra>.faa
```

Quando presentes, também aproveita:

```text
<amostra>.inference.tsv
<amostra>.hypotheticals.faa
<amostra>.hypotheticals.tsv
```

O JSON fornece metadados de anotação e o FAA é tratado como fonte contratual das sequências proteicas.

Antes de classificar as proteínas, o script verifica consistência de identificadores e sequências entre os arquivos.

---

## 18. Executar primeiro em modo `--dry-run`

Antes de criar qualquer saída de seleção:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --dry-run
```

O `--dry-run`:

- localiza as amostras;
- identifica o JSON e o FAA principais;
- informa os arquivos auxiliares encontrados;
- verifica os diretórios que seriam utilizados;
- não produz os arquivos finais de seleção.

Esse deve ser o primeiro comando utilizado nessa etapa.

---

## 19. Executar a seleção de todas as amostras

Depois de um `--dry-run` sem problemas:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT"
```

O programa percorre as amostras válidas e informa, para cada uma, as quantidades classificadas como:

```text
high
medium
low
sorf
pseudogene
```

Ao final, também informa os totais de amostras concluídas, ignoradas ou com falha.

---

## 20. Processar apenas uma ou algumas amostras

### Uma amostra

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --sample amostra_01
```

### Várias amostras

A opção `--sample` pode ser repetida:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --sample amostra_01 \
  --sample amostra_02
```

### Lista em arquivo

Também é possível fornecer um arquivo contendo uma amostra por linha:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --samples-file amostras.txt
```

Linhas iniciadas por `#` podem ser utilizadas como comentários na lista.

---

## 21. Modo direto para uma única pasta Bakta

Quando se deseja analisar diretamente uma pasta específica:

```bash
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-dir /caminho/output00_bakta/amostra_01 \
  --outdir /caminho/output00_bakta_selecao/amostra_01
```

Se o prefixo dos arquivos não puder ser inferido automaticamente, pode-se informar:

```bash
--prefix PREFIXO
```

---

## 22. Política padrão de classificação

Os valores padrão são:

```text
proteína curta:           < 90 aminoácidos
identidade fraca:         < 90%
cobertura da query fraca: < 80%
cobertura do alvo fraca:  < 80%
```

Esses limiares podem ser modificados com:

```text
--short-aa-threshold
--weak-identity
--weak-query-cov
--weak-subject-cov
```

Os três últimos recebem valores entre `0` e `1`.

Por exemplo, `0.90` corresponde a 90%.

### Prioridade `high`

Em geral, reúne proteínas com descrição ausente, hipotética ou não caracterizada.

### Prioridade `medium`

Em geral, reúne descrições genéricas em nível de família, domínio, DUF ou UPF.

### Prioridade `low`

Em geral, reúne descrições incertas ou anotações informativas sustentadas por identidade ou cobertura abaixo dos limiares definidos.

### Não selecionadas

Proteínas com produto informativo e evidência considerada adequada permanecem registradas na auditoria, mas não são enviadas aos FASTAs de prioridade.

---

## 23. Flags de auditoria não são automaticamente erros

O script pode registrar flags como:

```text
short_protein
gene_missing
ec_missing
cog_missing
go_missing
inference_missing
edge_feature
truncated_...
```

Essas flags descrevem características que podem ser úteis para revisão.

Elas **não significam, isoladamente, que a anotação esteja errada** e não são necessariamente critérios de seleção.

---

## 24. sORFs e pseudogenes

Por padrão:

- sORFs são preservadas separadamente;
- pseudogenes são preservados separadamente;
- nenhum dos dois grupos é misturado automaticamente aos alvos convencionais.

Existem opções para alterar esse comportamento:

```text
--include-sorf-targets
--include-pseudogenes
```

Essas opções devem ser utilizadas conscientemente, pois modificam o conjunto de proteínas submetido à classificação convencional.

---

## 25. Opções que exigem revisão prévia

O seletor foi projetado para interromper a execução quando encontra determinadas inconsistências.

As seguintes opções existem para permitir continuar **somente depois de uma revisão consciente do problema**:

```text
--allow-id-mismatch
--allow-sequence-mismatch
--allow-hypothetical-mismatch
```

Elas significam, respectivamente:

- aceitar diferenças entre os loci do JSON e do FAA;
- aceitar divergências entre as sequências do JSON e do FAA;
- aceitar inconsistências envolvendo o subconjunto oficial `hypotheticals.faa`.

Não utilize essas opções como forma de simplesmente eliminar uma mensagem de erro.

Primeiro determine a causa da inconsistência.

---

## 26. Uso de `--force`

Os dois scripts possuem situações em que existe uma opção `--force`, mas ela deve ser interpretada de acordo com a etapa.

### No `fleet`

```text
work --force
```

passa `--force` ao Bakta e permite sobrescrever uma pasta de saída já existente.

### No seletor

```text
--force
```

permite substituir apenas as saídas conhecidas da etapa de seleção.

Em condições normais, o seletor possui comportamento seguro de reexecução: resultados válidos podem ser reconhecidos pelo manifesto e pelos fingerprints dos arquivos de entrada e, nesse caso, são ignorados em vez de recriados.

---

## 27. Interromper o lote no primeiro erro

Por padrão, o seletor pode continuar para as outras amostras depois de uma falha individual.

Para interromper imediatamente:

```bash
--fail-fast
```

Isso é particularmente útil durante validação ou desenvolvimento.

---

## 28. Arquivos produzidos pela seleção

Para cada amostra, a etapa pode produzir:

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

### `targets_refined.tsv`

Contém as proteínas selecionadas, as informações originais do Bakta, os dados de inferência, as flags e a justificativa da prioridade.

### `selection_audit.tsv`

Registra a decisão tomada para as proteínas auditadas, inclusive proteínas não selecionadas ou tratadas separadamente.

### FASTAs de prioridade

```text
priority_high.faa
priority_medium.faa
priority_low.faa
```

contêm as sequências destinadas às buscas posteriores de acordo com a prioridade definida.

### `target_report.tsv`

Resume contagens gerais e frequência das principais flags.

### `target_manifest.json`

Registra parâmetros, versões, arquivos de entrada, hashes SHA-256, contagens, avisos e fingerprints das saídas.

É o principal registro de rastreabilidade da etapa de seleção.

---

## 29. Verificar os testes do repositório

Os testes básicos podem ser executados localmente sem rodar uma anotação Bakta completa.

Primeiro, verifique a sintaxe dos scripts:

```bash
python -m py_compile \
  scripts/anot00_bakta_01_fleet.py \
  scripts/anot00_bakta_02_selecionados.py
```

Depois execute os testes:

```bash
python -m unittest discover -s tests -v
```

Os testes verificam, entre outros pontos:

- reconhecimento do FASTA genômico principal;
- rejeição dos FASTAs de CDS e RNA;
- saneamento dos nomes das amostras;
- normalização de accession;
- classificação de exemplos `high`, `medium` e `low`;
- não seleção de uma anotação informativa com evidência forte.

O GitHub Actions executa automaticamente essas verificações em Python 3.9, 3.10, 3.11, 3.12 e 3.13 a cada atualização do repositório.

---

## 30. Sequência mínima recomendada

Para quem deseja apenas seguir o fluxo normal, a sequência é:

```bash
# 1. Entrar no repositório
cd bakta-distributed-annotation

# 2. Ativar o ambiente
mamba activate bakta_env

# 3. Criar a fila
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  init

# 4. Testar um job
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work --id worker-01 --cpus 4 --db "$BAKTA_DB" --skip-plot --once

# 5. Processar o restante da fila
python3 scripts/anot00_bakta_01_fleet.py \
  --input-root "$INPUT_ROOT" \
  --output-root "$BAKTA_OUTPUT" \
  work --id worker-01 --cpus 4 --db "$BAKTA_DB" --skip-plot --until-empty

# 6. Conferir a seleção sem escrever resultados
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT" \
  --dry-run

# 7. Executar a auditoria e seleção
python3 scripts/anot00_bakta_02_selecionados.py \
  --bakta-root "$BAKTA_OUTPUT" \
  --output-root "$SELECTION_OUTPUT"
```

---

## 31. Fluxo conceitual

```text
Pastas com genomas do NCBI
            ↓
     criação da fila
            ↓
 pending → running → done/failed
            ↓
 execução do Bakta por um ou mais workers
            ↓
       output00_bakta
            ↓
     dry-run da seleção
            ↓
 auditoria JSON × FAA × arquivos auxiliares
            ↓
 classificação high / medium / low
            ↓
 FASTAs + tabelas + relatório + manifesto
            ↓
 refinamento funcional posterior
```

---

## 32. Escopo do repositório

Este repositório termina na anotação com Bakta e na seleção/auditoria imediatamente posterior.

O arquivo:

```text
<amostra>.priority_high.faa
```

foi concebido como contrato de entrada para a próxima etapa curada do pipeline mais amplo, mas as etapas posteriores não fazem parte deste repositório.

Para detalhes sobre a finalidade dos scripts, consulte também o [`README.md`](../README.md).
