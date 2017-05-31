#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
# library(caret)

theme_set(theme_minimal(20))
essential_genes <- scan("/ifs/work/taylorlab/pensona/negative_selection/essential_genes.csv", what = "")
isoform_overrides_at_mskcc <- fread("/ifs/work/taylorlab/pensona/negative_selection/isoform_overrides_at_mskcc")
OncoKB_definitions <- fread("/ifs/work/taylorlab/pensona/negative_selection/OncoKB_definitions.txt")

genes_table <- as.data.table(readxl::read_excel("/ifs/work/taylorlab/pensona/negative_selection/nature12213-s2.xls", 5))
genes_table[, replication_time := as.integer(replication_time)]
genes_table[, log_noncoding_mutation_rate := log10(noncoding_mutation_rate)]
genes_table[, log_expression_CCLE := log10(expression_CCLE)]

my_fn <- function(data, mapping, method="gam", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(...)
  p
}
GGally::ggpairs(genes_table[, .(log_expression_CCLE, replication_time, HiC_compartment, local_GC_content, log_noncoding_mutation_rate)],
                lower = list(continuous = my_fn))

genes_table_filt <- genes_table[!is.na(replication_time) & !is.na(local_GC_content) & !is.na(HiC_compartment) & !is.infinite(log_expression_CCLE)]
x <- as.matrix(genes_table_filt[, .(log_expression_CCLE, replication_time, HiC_compartment, local_GC_content)])
y <- genes_table_filt$log_noncoding_mutation_rate
# fit <- glm.fit(x, y)
# fit <- glmnet::glmnet(x, y, family="gaussian")
set.seed(100)
fit <- caret::train(x, y, method = "glmnet", trControl = trainControl(method = "repeatedcv", number = 10, savePredictions = "final"))
pred <- setDT(fit$pred)[order(rowIndex)]
pred[, nmut_pred := 10^pred]
pred[, nmut_obs := 10^obs]
pred[, gene_name := genes_table_filt$gene]
pred[, essential := gene_name %in% essential_genes]
pred <- pred[order(essential)]
ggplot(pred, aes(nmut_pred, nmut_obs, col = essential)) + geom_text(aes(label = gene_name)) + scale_x_log10() + scale_y_log10() + annotation_logticks() + geom_abline(col = "black")

# ggplot(genes_table, aes(local_GC_content, noncoding_mutation_rate)) + geom_point() + scale_y_log10()
# ggplot(genes_table, aes(local_GC_content, noncoding_mutation_rate)) + geom_point() + scale_y_log10()


data_mutations_unfiltered <- fread("/ifs/work/taylorlab/pensona/molecular_diagnosis/msk-impact/msk-impact-clone/msk-impact/data_mutations_unfiltered.txt")
data_mutations_unfiltered[, nsyn := ifelse(Consequence %like% paste(autospy::Nonsyn_Consequences, collapse = "|"), "nsyn", ifelse(Consequence %like% "synonymous_variant", "syn", "other"))]
IMPACT341 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv3/genelist", what = "")
nmut_syn <- data_mutations_unfiltered[Hugo_Symbol %in% IMPACT341 & !Consequence %like% paste(autospy::Nonsyn_Consequences, collapse = "|") & Consequence %like% "synonymous_variant", .N, Hugo_Symbol]

translations_fasta <- Biostrings::readAAStringSet("/ifs/work/taylorlab/pensona/negative_selection/gencode.v19.pc_translations.fa.gz")
translations <- as.data.table(tstrsplit(names(translations_fasta), "\\|"))
setnames(translations,
         c("transcript_id", "gene_id",
           "havana_gene", "havana_transcript",
           "transcript_name", "gene_name",
           "transcript_length_aa"))
translations[, seq := as.character(translations_fasta)]
translations[, transcript_length_aa := as.integer(transcript_length_aa)]
translations[, transcript_id_stem := gsub("\\.[0-9]+$", "", transcript_id)]
# translations <- structure(as.character(proteins_fasta), .Names = proteins$protein_id)

isoform_overrides_at_mskcc <- fread("/ifs/work/taylorlab/pensona/negative_selection/isoform_overrides_at_mskcc")
isoform_overrides_at_mskcc <- isoform_overrides_at_mskcc[, 1:2, with = F]
setnames(isoform_overrides_at_mskcc, c("transcript_id_stem", "Hugo_Symbol"))
canonical_transcripts <- fread("/ifs/res/kandoth/dmp_vs_vep_canonical/vep_canonical", header = F)
canonical_transcripts <- canonical_transcripts[, 1:2, with = F]
setnames(canonical_transcripts, c("Hugo_Symbol", "transcript_id_stem"))
canonical_transcripts <- rbind(isoform_overrides_at_mskcc, canonical_transcripts)
canonical_transcripts <- canonical_transcripts[!duplicated(canonical_transcripts, by = "Hugo_Symbol")]
canonical_transcripts <- merge(canonical_transcripts, proteins, by = "transcript_id_stem")

nmut_syn <- merge(nmut_syn, canonical_transcripts[, .(Hugo_Symbol, transcript_length_bp)], by = "Hugo_Symbol")
nmut_syn[, Nrate := N / transcript_length_bp]

null_muts <- canonical_transcripts[gene_name %in% IMPACT341, .(Amino_Acid_Position = 1:transcript_length_aa, Nmut = 0), gene_name]

data_mutations_unfiltered[Variant_Type == "SNP", Amino_Acid_Position := as.integer(stringr::str_extract(HGVSp_Short, "(?<=^p.[A-Z])[0-9]+"))]


count_Amino_Acid_Position <- function(gene_name = "KRAS", 
                                      data_mutations_unfiltered = data_mutations_unfiltered, 
                                      nsyn = "nsyn",
                                      N_AA = 189,
                                      N_samples = 17816,
                                      N_bkg = 0.5){
  muts <- merge(
    data_mutations_unfiltered[
      Hugo_Symbol == gene_name & !is.na(Amino_Acid_Position), 
      .(Hugo_Symbol = unique(Hugo_Symbol), 
        N_mut = uniqueN(Tumor_Sample_Barcode)), 
      by = Amino_Acid_Position], 
    data.table(Hugo_Symbol = gene_name, 
               Amino_Acid_Position = 1:N_AA),
    all.y = T,
    by = c("Hugo_Symbol", "Amino_Acid_Position"))  
  muts[is.na(N_mut), N_mut := 0]
  
  set.seed(100)
  muts[, Amino_Acid_Rank := frank(-N_mut, ties.method = "random")]
  muts <- muts[order(Amino_Acid_Rank)]
  
  copy(muts)
}

plot_rank <- function(muts){
  if(!"expt" %in% names(muts)) muts[, expt := "obs"]
  set.seed(100)
  muts[, Amino_Acid_Rank := frank(-N_mut, ties.method = "random"), expt]
  muts <- muts[order(Amino_Acid_Rank)]
  
  muts[, N_samples := N_samples]
  muts[, paste0("mutrate", c("", "_lower", "_upper")) := as.list(
    binom::binom.confint(N_mut, N_samples, method = "exact",
                         conf.level = pnorm(1) - pnorm(-1))[4:6]), 
    .(N_mut, N_samples)]
  
  ggplot() +
    geom_point(data = muts[expt != "obs"], 
               aes(Amino_Acid_Rank, 
                   y = mutrate, 
                   ymin = mutrate_lower,
                   ymax = mutrate_upper,
                   col =  expt),
               alpha = 0.5) + 
    geom_point(data = muts[expt == "obs"], 
               aes(Amino_Acid_Rank, 
                   y = mutrate, 
                   ymin = mutrate_lower,
                   ymax = mutrate_upper),
               col =  "black",
               alpha = 1) + 
    #geom_errorbar(alpha = 0.3) + 
    scale_x_log10() + 
    scale_y_log10() + 
    facet_wrap(~Hugo_Symbol)
}


fit_counts <- function(muts,
                       N_AA = 189,
                       N_samples = 17816,
                       N_bkg = 0.5,
                       N_rm_obs = 7,
                       Hugo_Symbol = "TP53"){
  global_N_samples <- N_samples
  muts[, N_samples := global_N_samples]
  bbfits <- sapply(USE.NAMES = T, simplify = F,
    paste(N_rm_obs), 
    function(N_rm_obs){
      N_rm_obs <- as.integer(N_rm_obs)
      weights <- c(rep(1.0e-8, N_rm_obs),
                   rep(1, N_AA - N_rm_obs))
      
      set.seed(100)
      fit <- VGAM::vglm(cbind(N_mut, N_samples-N_mut) ~ 1, 
                        betabinomial, 
                        data = muts, 
                        weights = weights,
                        trace = TRUE)
      print(Coef(fit))
      print(fit@criterion$loglikelihood)
      
      fit_res <- as.data.table(transpose(
        lapply(
          transpose(
            replicate(100, 
                      simplify = F, 
                      sort(VGAM::rbetabinom(
                        N_AA, 
                        N_samples, 
                        Coef(fit)[[1]], 
                        Coef(fit)[[2]])))), 
          function(x) c(mean = mean(x), sd = sd(x)))))
      setnames(fit_res, c("N_mut", "sd_mut"))
      fit_res[, Amino_Acid_Rank := N_AA:1]
      fit_res[, Hugo_Symbol := Hugo_Symbol]
      fit_res[, mu := Coef(fit)[[1]]]
      fit_res[, rho := Coef(fit)[[2]]]
      fit_res[, loglikelihood := fit@criterion$loglikelihood]
      fit_res
    })
  
  # bbfits = sapply(USE.NAMES = T, simplify = F,
  #                 fits, 
  #                 function(fit){
  #                   data.table(
  #                     Hugo_Symbol = gene_name, 
  #                     Amino_Acid_Position = 1:N_AA, 
  #                     c("N_mut", "sd_mut") = transpose(
  #                       lapply(
  #                         transpose(
  #                           replicate(10, 
  #                                     simplify = F, 
  #                                     sort(VGAM::rbetabinom(
  #                                       N_AA, 
  #                                       N_samples, 
  #                                       Coef(fit)[[1]], 
  #                                       Coef(fit)[[2]])))), 
  #                         function(x) c(mean = mean(x), sd = sd(x)))))
  #                 })
  muts <- rbindlist(
    idcol = "expt",
    fill = TRUE,
    c(bbfits, list(obs = muts))
  )
  muts[, N_samples := global_N_samples]
  copy(muts)
}

fit_counts_ZIBB <- function(muts,
                            N_AA = 189,
                            N_samples = 17816,
                            N_bkg = 0.5,
                            N_rm_obs = 7,
                            Hugo_Symbol = "TP53"){
  global_N_samples <- N_samples
  muts[, N_samples := global_N_samples]
  bbfits <- sapply(USE.NAMES = T, simplify = F,
                   paste(N_rm_obs), 
                   function(N_rm_obs){
                     N_rm_obs <- as.integer(N_rm_obs)
                     weights <- c(rep(0, N_rm_obs),
                                  rep(1, N_AA - N_rm_obs))
                     
                     set.seed(100)
                     fit <- gamlss(cbind(N_mut, N_samples - N_mut) ~ 1, 
                                   family=ZIBB, 
                                   data=muts, 
                                   weights = weights,
                                   method = RS(), 
                                   control = gamlss.control(
                                     n.cyc = 100,
                                     trace = F))
                     params <- sapply(predictAll(fit), function(x)x[[1]])[c("mu", "sigma", "nu")]
                     fit_res <- as.data.table(transpose(
                       lapply(
                         transpose(
                           replicate(3, 
                                     simplify = F, 
                                     sort(do.call("rZIBB", c(n=N_AA, bd = N_samples, as.list(params)))))), 
                         function(x) c(mean = mean(x), sd = sd(x)))))
                     setnames(fit_res, c("N_mut", "sd_mut"))
                     fit_res[, Amino_Acid_Rank := N_AA:1]
                     fit_res[, Hugo_Symbol := Hugo_Symbol]
                     # fit_res[, mu := Coef(fit)[[1]]]
                     # fit_res[, rho := Coef(fit)[[2]]]
                     # fit_res[, loglikelihood := fit@criterion$loglikelihood]
                     fit_res
                   })
  
  # bbfits = sapply(USE.NAMES = T, simplify = F,
  #                 fits, 
  #                 function(fit){
  #                   data.table(
  #                     Hugo_Symbol = gene_name, 
  #                     Amino_Acid_Position = 1:N_AA, 
  #                     c("N_mut", "sd_mut") = transpose(
  #                       lapply(
  #                         transpose(
  #                           replicate(10, 
  #                                     simplify = F, 
  #                                     sort(VGAM::rbetabinom(
  #                                       N_AA, 
  #                                       N_samples, 
  #                                       Coef(fit)[[1]], 
  #                                       Coef(fit)[[2]])))), 
  #                         function(x) c(mean = mean(x), sd = sd(x)))))
  #                 })
  muts <- rbindlist(
    idcol = "expt",
    fill = TRUE,
    c(bbfits, list(obs = muts))
  )
  muts[, N_samples := global_N_samples]
  copy(muts)
}

plot_hist <- function(muts, max_mut = NA){
  if(!"expt" %in% names(muts)) muts[, expt := "obs"]
  mut_counts <- muts[, .(N_AA_mut = .N), keyby = .(expt, N_mut = round(N_mut), Hugo_Symbol)]
  mut_counts[, N_AA := N_AA]
  mut_counts[, paste0("f_AA_mut", c("", "_lower", "_upper")) := as.list(
    binom::binom.confint(N_AA_mut, N_AA, method = "exact",
                         conf.level = pnorm(1) - pnorm(-1))[4:6]), 
    .(N_AA_mut, N_AA)]
  
  
  if(is.na(max_mut)) max_mut <- N_AA
  p <- ggplot()
  if(mut_counts[expt != "obs", .N] > 0){
    p <- p + geom_step(
      data = mut_counts[expt != "obs"],
      aes(x = N_mut - 0.5,
          y = f_AA_mut,
          col = expt),
      size = 1)
  }
  p <- p + 
    geom_errorbar(data = mut_counts[expt == "obs"],
                  aes(x = N_mut,
                      ymin = f_AA_mut_lower,
                      ymax = f_AA_mut_upper),
                  size = 1) +
    geom_point(data = mut_counts[expt == "obs"], 
               aes(x = N_mut, 
                   y = f_AA_mut), 
               size = 3) +
    # scale_x_continuous(breaks = 0:10*2) +
    scale_y_log10() + 
    coord_cartesian(xlim = c(-0.5, max_mut + 0.5)) +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~Hugo_Symbol)
  p
}
