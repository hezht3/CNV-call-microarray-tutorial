###########################################################################
# Compare CNV calling results from PennCNV and QuantiSNP
# Author: Johnathan He
# Date: 7/21/2022
###########################################################################

setwd("/users/zhe/GEARS_CNV/")
require(tidyverse)
require(ggvenn)

###########################################################################
# Data preprocessing
###########################################################################
# Import CNV calling results from PennCNV & QuantiSNP (GC model adjusted)
cnv_penn <- read_table("./PennCNV-1.0.5/OUTPUT/sampleall.adjusted.rawcnv",
                       col_names = c("region", "numsnp", "length", "cn", "sample", "start_snp", "end_snp"))

cnv_quan <- list.files("./quantisnp/OUTPUT/", pattern = "*.cnv") %>% 
    map_dfr(~ read_table(paste("./quantisnp/OUTPUT/", .x, sep = ""), skip = 1,
                         col_names = c("sample", "chr", "start_pos", "end_pos", "start_snp", "end_snp",
                                       "length", "probes", "cn",
                                       "BF_max", "BF_1", "BF_2", "BF_3", "BF_4", "BF_5", "BF_6"))) %>% 
    filter(probes > 1)

# Harmonize data stratucture
cnv_penn <- cnv_penn %>% 
    separate(region, c("chr", "start_pos", "end_pos")) %>% 
    mutate(chr = str_remove_all(chr, pattern = "chr")) %>% 
    mutate(numsnp = str_remove_all(numsnp, pattern = "numsnp=")) %>% 
    mutate(length = str_remove_all(length, pattern = "length=")) %>% 
    mutate(cn = substr(cn, 11, 11)) %>% 
    mutate(cn = as.numeric(cn)) %>% 
    mutate(start_snp = str_remove_all(start_snp, pattern = "startsnp=")) %>% 
    mutate(end_snp = str_remove_all(end_snp, pattern = "endsnp=")) %>% 
    mutate(sample = case_when(sample == "INPUT/sample1.txt" ~ "99HI0697A",
                              sample == "INPUT/sample2.txt" ~ "99HI0698C",
                              sample == "INPUT/sample3.txt" ~ "99HI0700A"))

cnv_quan <- cnv_quan %>% 
    mutate(across(c("chr", "start_pos", "end_pos"),
                  ~ as.character(.x)))

# Combine calls
cnv_combine <- cnv_penn %>% 
    full_join(cnv_quan, by = c("sample", "chr", "start_snp", "end_snp"),
              suffix = c(".penn", ".quan")) %>% 
    arrange(sample, chr, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan) %>% 
    select(sample, chr, numsnp, probes, start_snp, end_snp, cn.penn, cn.quan,
           length.penn, length.quan, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan, contains("BF")) %>% 
    rename("probes.penn" = "numsnp", "probes.quan" = "probes")

###########################################################################
# Raw results comparison
###########################################################################
# Venn diagram
cnv_venn <- cnv_combine %>% 
    mutate(penn = if_else(is.na(length.penn) == FALSE, 1, 0),
           quan = if_else(is.na(length.quan) == FALSE, 1, 0)) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".")

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))

# Copy number
cnv_combine %>% 
    filter(is.na(cn.penn) == FALSE & is.na(cn.quan) == FALSE) %>% 
    mutate(across(c(cn.penn, cn.quan), ~ as.numeric(.x))) %>% 
    ggplot(aes(x = cn.penn, y = cn.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of copy number status for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# No. of SNPs/probes concordance
cnv_combine %>% 
    filter(is.na(probes.penn) == FALSE & is.na(probes.quan) == FALSE) %>% 
    mutate(across(contains("probes"), ~ as.numeric(.x))) %>% 
    ggplot(aes(x = probes.penn, y = probes.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of no. of SNPs/probes for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Length concordance
cnv_combine %>% 
    filter(is.na(length.penn) == FALSE & is.na(length.quan) == FALSE) %>% 
    mutate(across(contains("length"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    ggplot(aes(x = length.penn, y = length.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of length for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Start position concordance
cnv_combine %>% 
    filter(is.na(start_pos.penn) == FALSE & is.na(start_pos.quan) == FALSE) %>% 
    mutate(across(contains("start_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    ggplot(aes(x = start_pos.penn, y = start_pos.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of start position for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# End position concordance
cnv_combine %>% 
    filter(is.na(end_pos.penn) == FALSE & is.na(end_pos.quan) == FALSE) %>% 
    mutate(across(contains("end_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    ggplot(aes(x = end_pos.penn, y = end_pos.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of end position for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

###########################################################################
# Not harmonized calls (same start_snp, different end_snp | same end_snp, different start_snp)
###########################################################################

cnv_combine %>% 
    group_by(sample, start_snp) %>% 
    filter(n() > 1)   # 7 CNV calls

cnv_combine %>% 
    group_by(sample, end_snp) %>% 
    filter(n() > 1)   # 5 CNV calls

cnv_aggregate <- cnv_combine %>% 
    mutate(penn = if_else(is.na(length.penn) == FALSE, 1, 0),
           quan = if_else(is.na(length.quan) == FALSE, 1, 0)) %>%
    group_by(sample, start_snp) %>% 
    mutate(across(c(penn, quan), ~ ifelse(n() != 1, 1, .x))) %>% 
    slice(1) %>% 
    ungroup() %>% 
    group_by(sample, end_snp) %>% 
    mutate(across(c(penn, quan), ~ ifelse(n() != 1, 1, .x))) %>% 
    slice(1) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".")   # arbitarily keep 1 record per not harmonized CNV calls

ggvenn(list(`PennCNV` = cnv_aggregate %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_aggregate %>% filter(quan == 1) %>% pull(cnv)))

###########################################################################
# Consecutive probes > 5
###########################################################################

cnv_aggregate_qc <- cnv_aggregate %>% 
    filter((probes.penn > 5 & is.na(probes.quan) == TRUE) | (is.na(probes.penn) == TRUE & probes.quan > 5) |
               (probes.penn > 5 & probes.quan > 5))

ggvenn(list(`PennCNV` = cnv_aggregate_qc %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_aggregate_qc %>% filter(quan == 1) %>% pull(cnv)))
