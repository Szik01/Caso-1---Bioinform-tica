# Caso-1---Bioinformática

# Paciente com síndrome respiratória de origem desconhecida
Descrição: Paciente brasileira, sexo feminino, 21 anos, atendida no pronto atendimento de um
grande hospital em São Paulo com febre alta, tosse seca, dispneia e hipoxemia. Tomografia de
tórax evidenciando opacidades em vidro fosco bilaterais e histórico de viagem recente a Hong
Kong, na China. Foi feito o painel respiratório para múltiplos patógenos conhecidos, que veio
negativo (plot do filme Contágio).
Procedimento: Metagenômica de DNA total (shotgun) em amostra respiratória swab nasofaríngeo
solicitado pelo Ministério da Saúde.
Objetivo: Identificar agente etiológico causador da doença assim como emitir recomendações ao
Ministério da Saúde.

## Planejamento do Pipeline Metagenômico:
0) Coleta da amostra do paciente.
1) Preparação do ambiente padrão computacional.
2) Controle de qualidade.
3) Preparar o banco de dados para busca de patógenos.
4) Remoção de contaminantes.
5) Classificação Taxonômica.
6) Anotação e validação dos achados.
7) Laudamento.  

## Código Fonte e resultados (Realizados no Google Collab):

# 0) Configuração do ambiente padrão

%%bash
"Instalação miniconda e mamba"
wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /usr/local -u
mamba init
mamba config set auto_activate_base true
mamba shell init --shell bash --root-prefix=~/.local/share/mamba
rm Miniforge3-Linux-x86_64.sh

%%bash
"Instalação dos softwares necessários"
mamba install -c bioconda -c conda-forge fastqc bwa cutadapt samtools kraken2 krona blast seqtk spades -y

%%bash
"Criar as pastas para organizar os arquivos que serão gerados"
mkdir raw fastqc cutadapt bwa krakenDB kraken2 krona blast spades

%%bash
"Mover os arquivos para as respectivas pastas"
mv *.fastq.gz raw/
gunzip -c Chr15-reference.fasta.gz | tee bwa/Chr15-reference.fasta > blast/Chr15-reference.fasta

# 1) Controle de qualidade

%%bash
"Gerar o relatório da qualidade dos reads de sequenciamento"
fastqc -o fastqc/ raw/caso1-sindrome-respiratoriav2_R1.fastq.gz raw/caso1-sindrome-respiratoriav2_R2.fastq.gz

%%bash
"Filtrar sequências muito pequenas e cortar as pontas"
cutadapt \
  -u 5 -U 5 \
  -u -5 -U -5 \
  -m 90 \
  -o cutadapt/caso1-sindrome-respiratoriav2-trimmed_R1.fastq.gz \
  -p cutadapt/caso1-sindrome-respiratoriav2-trimmed_R2.fastq.gz \
  raw/caso1-sindrome-respiratoriav2_R1.fastq.gz \
  raw/caso1-sindrome-respiratoriav2_R2.fastq.gz

# 3) Preparar o banco de dados

%%bash
"Download banco de bactérias, fungos, vírus e humano"
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20250714.tar.gz

%%bash
"Extrair o banco para uma pasta"
tar -xzvf k2_pluspf_08_GB_20250714.tar.gz -C krakenDB/
rm k2_pluspf_08_GB_20250714.tar.gz

%%bash
"Download do banco de definições de taxonomia do NCBI"
bash /usr/local/opt/krona/updateTaxonomy.sh

# 4) Identificação Taxonômica

%%bash
"Buscar cada sequência contra um banco de dados e fazer a classificação"
kraken2 --paired cutadapt/caso1-sindrome-respiratoriav2-trimmed_R1.fastq.gz cutadapt/caso1-sindrome-respiratoriav2-trimmed_R2.fastq.gz \
-db krakenDB/ --report kraken2/caso1-sindrome-respiratoriav2-report.txt --report-minimizer --output kraken2/caso1-sindrome-respiratoriav2-output.txt

%%bash
"Gerar um gráfico de pizza iterativo com a classificação do kraken"
ktImportTaxonomy kraken2/caso1-sindrome-respiratoriav2-output.txt \
-o krona/caso1-sindrome-respiratoriav2-kraken.html -q 2 -t 3

# Com isso, é plotado os primeiros resultados do teste metagenômico:
Ele apresenta um gráfico de pizza interativo no qual mostra todos os tipos de células que compõe a coleta do exame, sendo elas iguais a 81% de organismos celulares, 7% de bactérias, 11% não foram identificados e por fim, 0,002% de vírus (Representado pelo primeiro plot). Destrinchando o gráfico, percebe-se que as bactérias são constituídas, em sua maioria, por bactérias normais que já residem no organismo do ser humano (Representado pelo segundo plot). Sendo assim, foram descartadas as hipóteses de patologias envolvendo bactérias. Por fim, os vírus foram observados como um potencial patógeno para a doença do caso, pois demonstram vírus como o SARS Coronavírus e BatCoronavírus (Representado pelo terceiro plot). Porém, para não tomar decisões precipitadas, é necessário verificar o que os 11% que não foram identificados significam. Portanto, será necessário retirar esse contaminante humano e células que ja conhecemos o genoma para focar na identificação das células não identificadas.

# Primeiro plot com todas células:
<img width="780" height="557" alt="image" src="https://github.com/user-attachments/assets/0a3e3c1f-49dc-49c2-a73b-377973d132c0" />

# Bactérias encontradas no primeiro plot:
<img width="893" height="559" alt="image" src="https://github.com/user-attachments/assets/e36ad0d9-a0b6-4fe6-af8e-ea63f0443fd6" />

# Vírus encontrados no primerio plot: 
<img width="896" height="625" alt="image" src="https://github.com/user-attachments/assets/93965bd8-13e0-4c90-b2c9-d3f423914f36" />

# 5) Remoção de contaminantes humanos

%%bash
"Indexar o genoma humano de referência para realizar o mapeamento e posterior remoção de contaminante"
bwa index bwa/Chr15-reference.fasta.gz

%%bash
"Proceder com o alinhamento contra o genoma humano"
bwa mem bwa/Chr15-reference.fasta.gz \
  cutadapt/caso1-sindrome-respiratoriav2-trimmed_R1.fastq.gz \
  cutadapt/caso1-sindrome-respiratoriav2-trimmed_R2.fastq.gz \
  | samtools view -b > bwa/caso1-sindrome-respiratoriav2-mapped-host.bam

%%bash
"Com a informação de mapeamento contida no BAM, vamos filtrar somente o que não mapeeou" 
samtools view -u -f 12 \
  -b bwa/caso1-sindrome-respiratoriav2-mapped-host.bam \
  | samtools sort -n > bwa/caso1-sindrome-respiratoriav2-unmapped-host.bam

%%bash
"Gerar FASTQ de reads não mapeados"
samtools fastq bwa/caso1-sindrome-respiratoriav2-unmapped-host.bam \
  -1 bwa/caso1-sindrome-respiratoriav2-unmapped-host_R1.fastq \
  -2 bwa/caso1-sindrome-respiratoriav2-unmapped-host_R2.fastq

# 6) Montagem dos reads em sequências menores 

%%bash
"Realiza a junção das sequências menores em contigs maiores (montagem)"
spades.py --meta -1 bwa/caso1-sindrome-respiratoriav2-unmapped-host_R1.fastq -2 bwa/caso1-sindrome-respiratoriav2-unmapped-host_R2.fastq -o spades/

# 7) Nova busca

%%bash
"Configura o banco novo para buscas via BLAST"
makeblastdb -in blast/Chr15-reference.fasta -dbtype nucl

%%bash
"Realizar nova busca, agora com BLAST"
blastn -query spades/scaffolds.fasta -db blast/Chr15-reference.fasta -out blast/caso1-sindrome-respiratoriav2-blast-scaffolds.txt -outfmt 6

%%bash
"Gera o gráfico de pizza para visualização dos resultados"
ktImportTaxonomy blast/caso1-sindrome-respiratoriav2-blast-scaffolds.txt -o krona/caso1-sindrome-respiratoriav2-blast.html

FIM DO CÓDIGO FONTE ( pode ser acessado por esse link: https://colab.research.google.com/drive/1YIlkcQefToWCGpwYJ5wBFNXyR6X7v1Uc?usp=sharing) 

# Resultado: 

O gráfico plotado após a nova busca não apresenta nenhum resultado, como podemos observar abaixo. 

<img width="1286" height="631" alt="image" src="https://github.com/user-attachments/assets/4559fdf5-d06d-4c8f-b31b-4fb74ba6d702" />

Com isso, podemos realizar algumas afirmações sobre o caso abordado. 
* As leituras remanescentes não apresentaram alinhamento significativo com genomas humanos, bacterianos ou virais catalogados nas bases de referência;
* A montagem metagenômica (SPAdes) resultou em um conjunto limitado de contigs, sem similaridade substancial em BLASTn contra bancos genômicos atualizados;
* A análise taxonômica pelo Krona não retornou categorias atribuíveis, reforçando a ausência de correspondência com organismos previamente conhecidos.

Tal achado é altamente sugestivo da presença de um agente genético ainda não catalogado, possivelmente um vírus emergente com identidade parcial, porém insuficiente, com as famílias de coronavírus já conhecidas (SARS-CoV, MERS-CoV, e coronavírus de morcegos).
O comportamento observado — desaparecimento de classificações no Krona após filtragem, aliado à inexistência de hits significativos em BLAST — é compatível com a detecção de um novo agente etiológico, cujas sequências apresentam divergência genômica relevante em relação aos vírus atualmente disponíveis nos bancos de dados públicos.
