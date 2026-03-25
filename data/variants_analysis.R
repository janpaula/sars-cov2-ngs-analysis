library(ggplot2)
library(knitr)
library(VariantAnnotation)
library(GenomicFeatures)
library(GenomicRanges)
library(Gviz)
library(Biostrings)
library(txdbmaker)
library(rtracklayer)
library(VariantAnnotation)
library(ggplot2)

# Ler as mutações
mut <- read.table("mutation_positions.tsv")
colnames(mut) <- c("pos")

# Definir genes do genoma viral
genes <- data.frame(
  gene=c("ORF1ab","Spike","ORF3a","E","M","N"),
  start=c(266,21563,25393,26245,26523,28274),
  end=c(21555,25384,26220,26472,27191,29533)
)

# Classificar mutações no Spike

mut$region <- "Other"
mut$region[
  mut$pos >= 21563 &
    mut$pos <= 25384
] <- "Spike"

# Gráfico estilo artigo
ggplot() +
  
  geom_rect(
    data=genes,
    aes(xmin=start, xmax=end, ymin=0, ymax=0.5),
    fill="lightblue",
    color="black"
  ) +
  
  geom_text(
    data=genes,
    aes(x=(start+end)/2, y=0.7, label=gene),
    size=3
  ) +
  
  geom_point(
    data=mut,
    aes(x=pos, y=1, color=region),
    size=3
  ) +
  
  scale_color_manual(values=c(
    "Spike"="red",
    "Other"="black"
  )) +
  
  labs(
    title="Mutations across SARS-CoV-2 genome",
    x="Genome position (bp)",
    y=""
  ) +
  
  theme_minimal()

# Total de reads
stats <- readLines("mapping_stats.txt")

total_reads <- as.numeric(strsplit(stats[1], " ")[[1]][1])
mapped_reads <- as.numeric(strsplit(stats[5], " ")[[1]][1])

mapping_rate <- round(mapped_reads/total_reads*100,2)


# Cobertura média
mean_cov <- as.numeric(
  gsub(",", ".", readLines("mean_coverage.txt"))
)

# Número de variantes
variant_count <- scan("variant_count.txt")

# Criar tabela resumo
summary_table <- data.frame(
  
  Metric = c(
    "Total reads",
    "Mapped reads",
    "Mapping rate (%)",
    "Mean genome coverage",
    "Number of variants"
  ),
  
  Value = c(
    total_reads,
    mapped_reads,
    mapping_rate,
    round(mean_cov,2),
    variant_count
  )
  
)


# Figura estilo Nature / Cell
library(ggplot2)

genes <- data.frame(
  gene=c("ORF1ab","Spike","ORF3a","E","M","N"),
  start=c(266,21563,25393,26245,26523,28274),
  end=c(21555,25384,26220,26472,27191,29533)
)

mut <- read.table("vcf_results/mutation_positions.tsv")
colnames(mut) <- "pos"

ggplot() +
  
  geom_rect(
    data=genes,
    aes(xmin=start,xmax=end,ymin=0,ymax=0.4),
    fill="skyblue"
  ) +
  
  geom_text(
    data=genes,
    aes(x=(start+end)/2,y=0.5,label=gene),
    size=3
  ) +
  
  geom_point(
    data=mut,
    aes(x=pos,y=1),
    color="red",
    size=3
  ) +
  
  labs(
    title="Mutation landscape across the SARS-CoV-2 genome",
    x="Genome position (bp)",
    y=""
  ) +
  
  theme_classic()

# ANÁLISE DAS VARIANTES POR QUALIDADE
# Importar o arquivo VCF
vcf <- import("vcf_results/variants.vcf.gz")

# Converter para DataFrame
vcf_df <- as.data.frame(rowRanges(vcf))
head(vcf_df)

# Filtrar variantes por qualidade, removendo variantes com baixa qualidade.
vcf_highqual <- subset(vcf_df, QUAL > 30)

# Filtro mantém apenas variantes com alta confiança estatística.
nrow(vcf_df)
nrow(vcf_highqual)

# Verificar distribuição de Qualidade
summary(vcf_df$QUAL)

# Visualizar variantes no genoma
ggplot(vcf_highqual,
       aes(x=start,
           y=QUAL)) +
  geom_point(size=3,color="red") +
  labs(
    title="High-confidence variants across the genome",
    x="Genome position",
    y="Quality score"
  ) +
  theme_classic()

# Plot 2
ggplot(vcf_df, aes(x=QUAL)) +
  geom_histogram(bins=20, fill="orange") +
  labs(
    title="Variant Quality Distribution",
    x="QUAL score",
    y="Number of Variants"
  )

# Visualizar posição das variantes no genoma
ggplot(vcf_highqual, aes(x=start, y=QUAL)) +
  geom_point(size=3, color="red") +
  labs(
    title="High-quality variants across the SARS-CoV-2 genome",
    x="Genome position",
    y="Variant quality"
  ) +
  theme_classic()

# Exportar tabela filtrada
write.table(vcf_highqual,
            "high_quality_variants.tsv",
            sep="\t",
            quote=FALSE,
            row.names=FALSE)

# Verificar quantas variantes realmente existem no VCF
vcf <- readVcf("vcf_results/variants.vcf.gz")
nrow(vcf)

# Verificar Profundidade
dp <- info(vcf)$DP
length(dp)
summary(dp)

# Verificar Frequencia Alélica
info(vcf)

af <- info(vcf)$AF
length(af)
summary(af)

length(dp)
length(af)

# Transformando em Lista
af_numeric <- unlist(info(vcf)$AF)

# Verificação se foi td ok
length(af_numeric)
summary(af_numeric)

#Verificando a classificação
class(af)

# Gráfico de Frequencia Alélica
ggplot(data.frame(AF=af_numeric), aes(x=AF)) +
  geom_histogram(bins=20, fill="steelblue") +
  labs(
    title="Allele Frequency Distribution",
    x="Allele Frequency",
    y="Number of Variants"
  )

# Profundidade da leitura READ DEPTH Distribution
dp <- info(vcf)$DP

# Conversão para lista numerica
dp_numeric <- unlist(info(vcf)$DP)

# Grafico Plotado
ggplot(data.frame(DP=dp_numeric), aes(x=DP)) +
  geom_histogram(bins=20, fill="darkgreen") +
  scale_x_log10() +
  labs(
    title="Read Depth Distribution",
    x="Read Depth (log scale)",
    y="Number of Variants"
  )

rowRanges(vcf)

# Cria tabela
vcf_df <- data.frame(
  CHROM = as.character(seqnames(rowRanges(vcf))),
  POS = start(rowRanges(vcf)),
  REF = as.character(ref(vcf)),
  ALT = as.character(unlist(alt(vcf))),
  QUAL = qual(vcf)
)

head(vcf_df)
vcf_highqual <- vcf_df[vcf_df$QUAL > 30, ]
nrow(vcf_highqual)

options(ucscChromosomeNames = FALSE)

# Gera gráficos
genome_axis <- GenomeAxisTrack()

variant_track <- AnnotationTrack(
  start = vcf_highqual$POS,
  end = vcf_highqual$POS,
  chromosome = "NC_045512",
  genome = "SARSCoV2",
  name = "Variants",
  col = "red",
  fill = "red"
)


gene_track <- AnnotationTrack(
  start = c(266,21563,25393,26245,27202,27394,27894,28274),
  end   = c(21555,25384,26220,26472,27387,27759,28259,29533),
  chromosome = "NC_045512",
  id = c("ORF1ab","Spike","ORF3a","E","M","ORF6","ORF7a","N"),
  genome = "SARSCoV2",
  name = "Genes",
  fill = "lightblue"
)

plotTracks(
  list(genome_axis, gene_track, variant_track),
  from = 1,
  to = 29903
)

# Gera tabela 
options(ucscChromosomeNames = FALSE)
vcf_df <- data.frame(
  CHROM = as.character(seqnames(rowRanges(vcf))),
  POS   = start(rowRanges(vcf)),
  REF   = as.character(ref(vcf)),
  ALT   = as.character(unlist(alt(vcf))),
  QUAL  = qual(vcf)
)

dp <- unlist(info(vcf)$DP)
af <- unlist(info(vcf)$AF)

vcf_df$DP <- dp
vcf_df$AF <- af

vcf_highqual <- vcf_df[vcf_df$QUAL > 30, ]

# Gera gráficos
gene_track <- AnnotationTrack(
  start = c(266,21563,25393,26245,27202,27394,27894,28274),
  end   = c(21555,25384,26220,26472,27387,27759,28259,29533),
  chromosome = "NC_045512",
  id = c("ORF1ab","Spike","ORF3a","E","M","ORF6","ORF7a","N"),
  genome = "SARSCoV2",
  name = "Genes",
  fill = "lightblue",
  col = "darkblue"
)

# Cria Tabela
variant_track <- AnnotationTrack(
  start = vcf_highqual$POS,
  end   = vcf_highqual$POS,
  chromosome = "NC_045512",
  genome = "SARSCoV2",
  name = "Variants",
  col = "red",
  fill = "red"
)

af_track <- DataTrack(
  start = vcf_df$POS,
  end   = vcf_df$POS,
  data  = vcf_df$AF,
  chromosome = "NC_045512",
  genome = "SARSCoV2",
  type = "p",
  name = "Allele Frequency",
  col = "purple"
)

depth_track <- DataTrack(
  start = vcf_df$POS,
  end   = vcf_df$POS,
  data  = vcf_df$DP,
  chromosome = "NC_045512",
  genome = "SARSCoV2",
  type = "histogram",
  name = "Read Depth",
  col.histogram = "darkgreen",
  fill.histogram = "darkgreen"
)

#Cria gráficos
axis_track <- GenomeAxisTrack()
vcf_df$variant_name <- paste0(vcf_df$POS, "_", vcf_df$REF, ">", vcf_df$ALT)
vcf_highqual <- vcf_df[vcf_df$QUAL > 30, ]
variant_track <- AnnotationTrack(
  start = vcf_highqual$POS,
  end   = vcf_highqual$POS,
  chromosome = "NC_045512",
  genome = "SARSCoV2",
  id = vcf_highqual$variant_name,
  name = "Variants",
  col = "red",
  fill = "red"
)

plotTracks(
  list(axis_track, gene_track, variant_track, depth_track, af_track),
  from = 1,
  to = 29903,
  main = "Genomic Distribution of Variants in SARS-CoV-2",
  showFeatureId = TRUE
)


plotTracks(
  list(axis_track, gene_track, variant_track, depth_track, af_track),
  from = 1,
  to = 29903,
  main = "Distribution of Variants Across the SARS-CoV-2 Genome"
)
