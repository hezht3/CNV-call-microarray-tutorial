###########################################################################
# Compare CNV calling results from PennCNV and QuantiSNP
# Author: Johnathan He
# Date: 7/21/2022
###########################################################################

setwd("/users/zhe/GEARS_CNV/")
require(tidyverse)
require(ggvenn)
require(gtsummary)

###########################################################################
# I. First combine concordant calls, then anneal
###########################################################################

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

# Harmonize data
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
                              sample == "INPUT/sample3.txt" ~ "99HI0700A")) %>% 
    mutate(software = "PennCNV")

cnv_quan <- cnv_quan %>% 
    mutate(across(c("chr", "start_pos", "end_pos"),
                  ~ as.character(.x))) %>% 
    mutate(software = "QuantiSNP")

###########################################################################
# Combine completely concordant calls (same start SNP, same end SNP)
###########################################################################
# Combine calls
cnv_combine <- cnv_penn %>% 
    full_join(cnv_quan, by = c("sample", "chr", "start_snp", "end_snp"),
              suffix = c(".penn", ".quan")) %>% 
    arrange(sample, chr, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan) %>% 
    select(sample, chr, numsnp, probes, start_snp, end_snp, cn.penn, cn.quan,
           length.penn, length.quan, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan, 
           software.penn, software.quan, contains("BF")) %>% 
    rename("probes.penn" = "numsnp", "probes.quan" = "probes") %>% 
    mutate(across(c(cn.penn, cn.quan), ~ as.numeric(.x))) %>% 
    mutate(across(contains("probes"), ~ as.numeric(.x))) %>% 
    mutate(across(contains("length"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    mutate(across(contains("start_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    mutate(across(contains("end_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x)))

# Venn diagram
cnv_venn <- cnv_combine %>% 
    mutate(penn = if_else(software.penn == "PennCNV", 1, 0),
           quan = if_else(software.quan == "QuantiSNP", 1, 0)) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".")

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))

# Validity assessment - copy number
cnv_combine %>% 
    filter(is.na(cn.penn) == FALSE & is.na(cn.quan) == FALSE) %>% 
    ggplot(aes(x = cn.penn, y = cn.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of copy number status for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - no. of SNPs/probes concordance
cnv_combine %>% 
    filter(is.na(probes.penn) == FALSE & is.na(probes.quan) == FALSE) %>% 
    ggplot(aes(x = probes.penn, y = probes.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of no. of SNPs/probes for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - length concordance
cnv_combine %>% 
    filter(is.na(length.penn) == FALSE & is.na(length.quan) == FALSE) %>% 
    ggplot(aes(x = length.penn, y = length.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of length for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - start position concordance
cnv_combine %>% 
    filter(is.na(start_pos.penn) == FALSE & is.na(start_pos.quan) == FALSE) %>% 
    ggplot(aes(x = start_pos.penn, y = start_pos.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of start position for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - end position concordance
cnv_combine %>% 
    filter(is.na(end_pos.penn) == FALSE & is.na(end_pos.quan) == FALSE) %>% 
    ggplot(aes(x = end_pos.penn, y = end_pos.quan)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    coord_equal() +
    ggtitle("Concordance of end position for each CNV called") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Filter completely concordant calls
cnv_concor_com <- cnv_combine %>%
    filter(if_all(contains("software"), ~ !is.na(.x))) %>% 
    select(sample, chr, start_snp, end_snp, contains("penn"), cn.quan, contains("BF"), - contains("software")) %>% 
    rename("probes" = "probes.penn", "length" = "length.penn", "start_pos" = "start_pos.penn", "end_pos" = "end_pos.penn")
    

###########################################################################
# Combine partially concordant calls (same start_snp, different end_snp | same end_snp, different start_snp)
###########################################################################
# Rules: generate unique cluster ID for partially concordant calls
# Step 1. Sort calls by sample ID, chromosome, starting position (1st), ending position (2nd)
# Step 2. Group calls by sample ID, chromosome, filter groups with >1 calls
# Step 3. Generate concordance indicator variable:
#    1. If starting position of current call <= cumulative maximum ending position of previous calls, indicator remains
#       same as previous row
#    2. If starting position of current call > cumulative maximum ending position of previous calls, indicator + 1
# Step 4. Generate CNV ID as sample ID & chromosome & cumulative sum of concordance indicator
# After this step, CNV ID is unique for each partially concordant call pair/group

# Illustration:
# +-----------------------+-----------------------------------+------------------------------+-----------------------+--------+
# |Rows                   |       CNV visualization           |       Logic                  | Concordance indicator | CNV ID |
# +-----------------------+-----------------------------------+------------------------------+-----------------------+--------+
# |1st row                | ==================                |                              | 1                     | 1      |
# |Later rows             |    ====                           | Start pos <= cummax(end pos) | 0                     | 1      |
# |                       |         =========                 | Start pos <= cummax(end pos) | 0                     | 1      |
# |                       |                    ======         | Start pos > cummax(end pos)  | 1                     | 2      |
# |                       |                      ======       | Start pos <= cummax(end pos) | 0                     | 2      |
# |                       |                             ===== | Start pos > cummax(end pos)  | 1                     | 3      |
# +-----------------------+-----------------------------------+------------------------------+-----------------------+--------+

# Filter partially concordant calls + singleton calls
cnv_concor_par <- cnv_combine %>% 
    filter(if_any(contains("software"), ~ is.na(.x))) %>% 
    # unite `probes`, `cn`, `length`, `start position`, `end position` outputed from both software
    unite("probes", contains("probes"), sep = "", na.rm = TRUE) %>% 
    unite("cn", contains("cn"), sep = "", na.rm = TRUE) %>% 
    unite("length", contains("length"), sep = "", na.rm = TRUE) %>% 
    unite("start_pos", contains("start_pos"), sep = "", na.rm = TRUE) %>% 
    unite("end_pos", contains("end_pos"), sep = "", na.rm = TRUE) %>% 
    unite("software", contains("software"), sep = "", na.rm = TRUE) %>% 
    mutate(across(c("chr", "probes", "cn", "length", "start_pos", "end_pos"), ~ as.numeric(.x))) %>% 
    group_by(sample, chr) %>% 
    arrange(sample, chr, start_pos, end_pos) %>% 
    mutate(across(contains("pos"), ~ as.numeric(.x)))

# Filter partially concordant calls
cnv_concor_par <- cnv_concor_par %>% 
    group_by(sample, chr) %>% 
    filter(n() > 1) %>% 
    mutate(id_num = 1:n()) %>% 
    mutate(cummax_end_pos = cummax(end_pos)) %>% 
    # already sorted by starting position
    # if current call not overlap with any previous call, `condi` + 1
    # if current call overlap with any previous call, `condi` remain the same
    mutate(con_indi = case_when(id_num == 1 ~ 1,
                                id_num > 1 & start_pos <= lag(cummax_end_pos) ~ 0,
                                id_num > 1 & start_pos > lag(cummax_end_pos) ~ 1)) %>% 
    mutate(con_indi_cumsum = cumsum(con_indi)) %>% 
    ungroup() %>% 
    mutate(cnv_id = paste(sample, chr, con_indi_cumsum, sep = ".")) %>%    # unique cluster id for partially concordant calls
    group_by(cnv_id) %>% 
    filter(n() > 1)
cnv_concor_par_cnvid <- cnv_concor_par %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% pull(cnv)

# Intersection
max2 <-  function(x) {
    u <- unique(x)
    sort(u, decreasing = TRUE)[2L]
}
min2 <-  function(x) {
    u <- unique(x)
    sort(u, decreasing = FALSE)[2L]
}
cnv_concor_par <- cnv_concor_par %>% 
    group_by(cnv_id) %>% 
    mutate(inter_start_pos = case_when(n() == 2 ~ max(start_pos),
                                       n() > 2 ~ max2(start_pos)),
           inter_end_pos = case_when(n() == 2 ~ min(end_pos),
                                     n() > 2 ~ min2(end_pos))) %>% 
    mutate(inter_length = inter_end_pos - inter_start_pos + 1) %>% 
    mutate(inter_length = if_else(n() == 3 & id_num == 2, lag(inter_length) + lead(inter_length), inter_length))

# Evaluate intersection regions - visualization
cnv_concor_par %>% distinct(cnv_id) %>% pull(cnv_id) %>% 
    map(~ cnv_concor_par %>% 
            filter(cnv_id == .x) %>% 
            mutate(cnv_call = paste(sample, chr, "call", id_num, software, sep = ".")) %>% 
            ggplot(aes(y = cnv_call, color = software)) +
            geom_linerange(aes(xmin = start_pos, xmax = end_pos), size = 3) +
            xlab("Position") +
            ylab("Sample ID.CHR.CNV.Software") +
            scale_color_manual("Software", values = c("#ff9900", "#146eb4")) +
            theme_minimal())

# Evaluate intersection regions - percent agreement [for 3 calls per CNV cluster, kept complete consesus region]
cnv_concor_eval <- cnv_concor_par$inter_length %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "sum") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "inter_length" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- inter_length),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE)

cnv_concor_eval %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    as.data.frame() %>% 
    select(perc.agree, software) %>% 
    mutate(perc.agree.60 = if_else(perc.agree >= 0.6, "Yes", "No")) %>% 
    mutate(perc.agree = as.character(round(perc.agree, 2))) %>% 
    mutate(software = factor(software,
                             levels = c("PennCNV", "QuantiSNP"),
                             labels = c("Gold standard: PennCNV", "Gold standard: QuantiSNP"))) %>% 
    tbl_summary(by = "software",
                label = list(perc.agree = "Percent Agreement",
                             perc.agree.60 = "Percent Agreement >= 60%")) %>% 
    bold_labels()

cnv_concor_eval %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    mutate(software = factor(software,
                             levels = c("PennCNV", "QuantiSNP"),
                             labels = c("Gold standard: PennCNV", "Gold standard: QuantiSNP"))) %>% 
    ggplot() +
    geom_histogram(aes(x = perc.agree, fill = software)) +
    facet_wrap(. ~ software, scales = "free_x") +
    geom_vline(aes(xintercept = 0.6, color = "red"), linetype = "dashed", show.legend = F) +
    xlab("Percent agreement") +
    ylab("Frequency") +
    scale_fill_manual("Software", values = c("#ff9900", "#146eb4"), guide = "none") +
    theme_minimal()

# Validity assessment - length concordance [for 3 calls per CNV cluster, kept complete consesus region]
cnv_concor_par$length %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "length" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- length),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = length.PennCNV, y = length.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of length for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(100, 210000) +
    ylim(100, 210000) +
    theme_minimal()

# Validity assessment - start position concordance [for 3 calls per CNV cluster, kept min start position]
cnv_concor_par$start_pos %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "min") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "start_pos" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- start_pos),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = start_pos.PennCNV, y = start_pos.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of start position for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - end position concordance [for 3 calls per CNV cluster, kept max start position]
cnv_concor_par$end_pos %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "end_pos" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- end_pos),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = end_pos.PennCNV, y = end_pos.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of end position for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - copy number [for 3 calls per CNV cluster, kept max cn relative to null 2]
cnv_concor_par$cn %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software),
              FUN = function(vector) {
                  if(min(vector) > 2) {return(max(vector))}
                  else if(max(vector) < 2) {return(min(vector))}
              }) %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "cn" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- cn),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = cn.PennCNV, y = cn.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of copy number status for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(0, 4) +
    ylim(0, 4) +
    theme_minimal()

# Validity assessment - no. of SNPs/probes concordance
cnv_concor_par$probes %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "probes" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- probes),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = probes.PennCNV, y = probes.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of no. of SNPs/probes for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(2, 11) +
    ylim(2, 11) +
    theme_minimal()

# Filter harmonized partially concordant calls [Pending: `NA` values filling in `probes`]
Mode <- function(x) {
    ux <- na.omit(unique(x) )
    tab <- tabulate(match(x, ux)); ux[tab == max(tab) ]
}
cnv_concor_par <- cnv_concor_par$cn %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software),
              FUN = function(vector) {
                  if(min(vector) > 2) {return(max(vector))}
                  else if(max(vector) < 2) {return(min(vector))}
              }) %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "cn" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- cn),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    group_by(cnv_id) %>% 
    mutate(inter_length = inter_end_pos - inter_start_pos + 1,
           inter_start_snp = ifelse(inter_start_pos == start_pos, start_snp, as.character(NA)),
           inter_end_snp = ifelse(inter_end_pos == end_pos, end_snp, as.character(NA)),
           inter_probes = ifelse(inter_start_pos == start_pos & inter_end_pos == end_pos, probes, as.numeric(NA))) %>% 
    mutate(across(contains("inter"), ~ Mode(.x))) %>% 
    mutate(across(contains("BF"), ~ Mode(.x))) %>% 
    as.data.frame() %>% 
    select(sample, chr, inter_start_snp, inter_end_snp, inter_probes, cn, inter_length, inter_start_pos, inter_end_pos, contains("BF"), software) %>% 
    pivot_wider(names_from = "software", values_from = "cn") %>% 
    rename_with(~ str_remove(., 'inter_')) %>% 
    rename("cn.penn" = "PennCNV", "cn.quan" = "QuantiSNP")

# Merge completely and partially concordant calls
cnv_concor <- cnv_concor_com %>% mutate(chr = as.numeric(chr)) %>% bind_rows(cnv_concor_par)

# Venn diagram
cnv_venn <- cnv_combine %>% 
    mutate(penn = if_else(software.penn == "PennCNV", 1, 0),
           quan = if_else(software.quan == "QuantiSNP", 1, 0)) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% 
    mutate(penn = if_else(cnv %in% cnv_concor_par_cnvid, 1, penn),
           quan = if_else(cnv %in% cnv_concor_par_cnvid, 1, quan))

cnv_venn <- cnv_venn %>% 
    filter(is.na(penn) | is.na(quan)) %>% 
    select(cnv, penn, quan) %>% 
    bind_rows(cnv_concor %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% mutate(penn = 1, quan = 1))

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))

# Export concordant calls for anneal
cnv_concor %>% 
    select(chr, start_pos, end_pos, cn.penn, sample, start_snp, end_snp, BF_max, probes) %>% 
    write.table("/users/zhe/GEARS_CNV/PennCNV-1.0.5/OUTPUT/sampleall.concordant.rawcnv",
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


###########################################################################
# II. First anneal, then combine concordant calls
###########################################################################

rm(list = ls())

# Evaluate anneal (PennCNV: 2, QuantiSNP: 0)
cnv_penn_pre <- read_table("./PennCNV-1.0.5/OUTPUT/sampleall.adjusted.rawcnv",
                           col_names = c("region", "numsnp", "length", "cn", "sample", "start_snp", "end_snp"))

cnv_penn_post <- read_table("./PennCNV-1.0.5/OUTPUT/sampleall.adjusted.anneal.rawcnv",
                            col_names = c("region", "numsnp", "length", "cn", "sample", "start_snp", "end_snp"))

cnv_penn <- list(cnv_penn_pre, cnv_penn_post) %>% 
    map(~ .x %>% 
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
                                      sample == "INPUT/sample3.txt" ~ "99HI0700A")) %>% 
            mutate(software = "PennCNV"))
cnv_penn_pre <- cnv_penn[[1]]
cnv_penn_post <- cnv_penn[[2]]

anneal_pre <- anti_join(cnv_penn_pre, cnv_penn_post, by = c("sample", "chr", "start_snp", "end_snp")) %>% 
    mutate(anneal = "pre")
anneal_post <- anti_join(cnv_penn_post, cnv_penn_pre, by = c("sample", "chr", "start_snp", "end_snp")) %>% 
    mutate(anneal = "post")

anneal_pre %>% 
    bind_rows(anneal_post) %>% 
    mutate(cnv_call = paste(sample, chr, start_pos, end_pos, sep = ".")) %>% 
    mutate(cnv_call = factor(cnv_call,
                             levels = c("99HI0698C.22.21361579.21392864", "99HI0698C.22.21401228.21554058", "99HI0698C.22.21361579.21554058"))) %>% 
    ggplot(aes(y = cnv_call, color = anneal)) +
    geom_linerange(aes(xmin = start_pos, xmax = end_pos), size = 3) +
    xlab("Position") +
    ylab("Sample ID.CHR.Start pos.End pos") +
    scale_color_manual("Anneal", values = c("#ff9900", "#146eb4")) +
    theme_minimal()

rm(list = ls())

###########################################################################
# Data preprocessing
###########################################################################
# Import annealed CNV calling results from PennCNV & QuantiSNP (GC model adjusted)

cnv_penn <- read_table("./PennCNV-1.0.5/OUTPUT/sampleall.adjusted.anneal.rawcnv",
                       col_names = c("region", "numsnp", "length", "cn", "sample", "start_snp", "end_snp"))
cnv_quan <- list.files("./quantisnp/OUTPUT/convert/", pattern = "*.adjust.penncnv") %>% 
    map_dfr(~ read_table(paste("./quantisnp/OUTPUT/convert/", .x, sep = ""),
                         col_names = c("region", "numsnp", "length", "cn", "sample", "start_snp", "end_snp")))

cnv_software <- list(cnv_penn, cnv_quan) %>% 
    map(~ .x %>% 
            separate(region, c("chr", "start_pos", "end_pos")) %>% 
            mutate(chr = str_remove_all(chr, pattern = "chr")) %>% 
            mutate(numsnp = str_remove_all(numsnp, pattern = "numsnp=")) %>% 
            mutate(length = str_remove_all(length, pattern = "length=")) %>% 
            mutate(cn = substr(cn, 11, 11)) %>% 
            mutate(cn = as.numeric(cn)) %>% 
            mutate(start_snp = str_remove_all(start_snp, pattern = "startsnp=")) %>% 
            mutate(end_snp = str_remove_all(end_snp, pattern = "endsnp=")))
cnv_penn <- cnv_software[[1]] %>%
    mutate(software = "PennCNV") %>% 
    mutate(sample = case_when(sample == "INPUT/sample1.txt" ~ "99HI0697A",
                              sample == "INPUT/sample2.txt" ~ "99HI0698C",
                              sample == "INPUT/sample3.txt" ~ "99HI0700A"))
cnv_quan <- cnv_software[[2]] %>%
    mutate(software = "QuantiSNP") %>%
    rename("probes" = "numsnp")

###########################################################################
# Combine completely concordant calls (same start SNP, same end SNP)
###########################################################################
# Combine calls
cnv_combine <- cnv_penn %>% 
    full_join(cnv_quan, by = c("sample", "chr", "start_snp", "end_snp"),
              suffix = c(".penn", ".quan")) %>% 
    arrange(sample, chr, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan) %>% 
    select(sample, chr, numsnp, probes, start_snp, end_snp, cn.penn, cn.quan,
           length.penn, length.quan, start_pos.penn, start_pos.quan, end_pos.penn, end_pos.quan, 
           software.penn, software.quan, contains("BF")) %>% 
    rename("probes.penn" = "numsnp", "probes.quan" = "probes") %>% 
    mutate(across(c(cn.penn, cn.quan), ~ as.numeric(.x))) %>% 
    mutate(across(contains("probes"), ~ as.numeric(.x))) %>% 
    mutate(across(contains("length"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    mutate(across(contains("start_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x))) %>% 
    mutate(across(contains("end_pos"), ~ str_remove_all(.x, pattern = ",") %>% as.numeric(.x)))

# Venn diagram
cnv_venn <- cnv_combine %>% 
    mutate(penn = if_else(software.penn == "PennCNV", 1, 0),
           quan = if_else(software.quan == "QuantiSNP", 1, 0)) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".")

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))

# Filter completely concordant calls
cnv_concor_com <- cnv_combine %>%
    filter(if_all(contains("software"), ~ !is.na(.x))) %>% 
    select(sample, chr, start_snp, end_snp, contains("penn"), cn.quan, contains("BF"), - contains("software")) %>% 
    rename("probes" = "probes.penn", "length" = "length.penn", "start_pos" = "start_pos.penn", "end_pos" = "end_pos.penn")

###########################################################################
# Combine partially concordant calls (same start_snp, different end_snp | same end_snp, different start_snp)
###########################################################################
# Filter partially concordant calls + singleton calls
cnv_concor_par <- cnv_combine %>% 
    filter(if_any(contains("software"), ~ is.na(.x))) %>% 
    # unite `probes`, `cn`, `length`, `start position`, `end position` outputed from both software
    unite("probes", contains("probes"), sep = "", na.rm = TRUE) %>% 
    unite("cn", contains("cn"), sep = "", na.rm = TRUE) %>% 
    unite("length", contains("length"), sep = "", na.rm = TRUE) %>% 
    unite("start_pos", contains("start_pos"), sep = "", na.rm = TRUE) %>% 
    unite("end_pos", contains("end_pos"), sep = "", na.rm = TRUE) %>% 
    unite("software", contains("software"), sep = "", na.rm = TRUE) %>% 
    mutate(across(c("chr", "probes", "cn", "length", "start_pos", "end_pos"), ~ as.numeric(.x))) %>% 
    group_by(sample, chr) %>% 
    arrange(sample, chr, start_pos, end_pos) %>% 
    mutate(across(contains("pos"), ~ as.numeric(.x)))

# Filter partially concordant calls
cnv_concor_par <- cnv_concor_par %>% 
    group_by(sample, chr) %>% 
    filter(n() > 1) %>% 
    mutate(id_num = 1:n()) %>% 
    mutate(cummax_end_pos = cummax(end_pos)) %>% 
    # already sorted by starting position
    # if current call not overlap with any previous call, `condi` + 1
    # if current call overlap with any previous call, `condi` remain the same
    mutate(con_indi = case_when(id_num == 1 ~ 1,
                                id_num > 1 & start_pos <= lag(cummax_end_pos) ~ 0,
                                id_num > 1 & start_pos > lag(cummax_end_pos) ~ 1)) %>% 
    mutate(con_indi_cumsum = cumsum(con_indi)) %>% 
    ungroup() %>% 
    mutate(cnv_id = paste(sample, chr, con_indi_cumsum, sep = ".")) %>%    # unique cluster id for partially concordant calls
    group_by(cnv_id) %>% 
    filter(n() > 1)
cnv_concor_par_cnvid <- cnv_concor_par %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% pull(cnv)

# Intersection
max2 <-  function(x) {
    u <- unique(x)
    sort(u, decreasing = TRUE)[2L]
}
min2 <-  function(x) {
    u <- unique(x)
    sort(u, decreasing = FALSE)[2L]
}
cnv_concor_par <- cnv_concor_par %>% 
    group_by(cnv_id) %>% 
    mutate(inter_start_pos = case_when(n() == 2 ~ max(start_pos),
                                       n() > 2 ~ max2(start_pos)),
           inter_end_pos = case_when(n() == 2 ~ min(end_pos),
                                     n() > 2 ~ min2(end_pos))) %>% 
    mutate(inter_length = inter_end_pos - inter_start_pos + 1) %>% 
    mutate(inter_length = if_else(n() == 3 & id_num == 2, lag(inter_length) + lead(inter_length), inter_length))

# Evaluate intersection regions - visualization
cnv_concor_par %>% distinct(cnv_id) %>% pull(cnv_id) %>% 
    map(~ cnv_concor_par %>% 
            filter(cnv_id == .x) %>% 
            mutate(cnv_call = paste(sample, chr, "call", id_num, software, sep = ".")) %>% 
            ggplot(aes(y = cnv_call, color = software)) +
            geom_linerange(aes(xmin = start_pos, xmax = end_pos), size = 3) +
            xlab("Position") +
            ylab("Sample ID.CHR.CNV.Software") +
            scale_color_manual("Software", values = c("#ff9900", "#146eb4")) +
            theme_minimal())

# Evaluate intersection regions - percent agreement [for 3 calls per CNV cluster, kept complete consesus region]
cnv_concor_eval <- cnv_concor_par$inter_length %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "sum") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "inter_length" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- inter_length),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE)

cnv_concor_eval %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    as.data.frame() %>% 
    select(perc.agree, software) %>% 
    mutate(perc.agree.60 = if_else(perc.agree >= 0.6, "Yes", "No")) %>% 
    mutate(perc.agree = as.character(round(perc.agree, 2))) %>% 
    mutate(software = factor(software,
                             levels = c("PennCNV", "QuantiSNP"),
                             labels = c("Gold standard: PennCNV", "Gold standard: QuantiSNP"))) %>% 
    tbl_summary(by = "software",
                label = list(perc.agree = "Percent Agreement",
                             perc.agree.60 = "Percent Agreement >= 60%")) %>% 
    bold_labels()

cnv_concor_eval %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    mutate(software = factor(software,
                             levels = c("PennCNV", "QuantiSNP"),
                             labels = c("Gold standard: PennCNV", "Gold standard: QuantiSNP"))) %>% 
    ggplot() +
    geom_histogram(aes(x = perc.agree, fill = software)) +
    facet_wrap(. ~ software, scales = "free_x") +
    geom_vline(aes(xintercept = 0.6, color = "red"), linetype = "dashed", show.legend = F) +
    xlab("Percent agreement") +
    ylab("Frequency") +
    scale_fill_manual("Software", values = c("#ff9900", "#146eb4"), guide = "none") +
    theme_minimal()

# Validity assessment - length concordance [for 3 calls per CNV cluster, kept complete consesus region]
cnv_concor_par$length %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "length" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- length),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = length.PennCNV, y = length.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of length for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(100, 210000) +
    ylim(100, 210000) +
    theme_minimal()

# Validity assessment - start position concordance [for 3 calls per CNV cluster, kept min start position]
cnv_concor_par$start_pos %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "min") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "start_pos" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- start_pos),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = start_pos.PennCNV, y = start_pos.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of start position for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - end position concordance [for 3 calls per CNV cluster, kept max start position]
cnv_concor_par$end_pos %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "end_pos" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- end_pos),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = end_pos.PennCNV, y = end_pos.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of end position for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    theme_minimal()

# Validity assessment - copy number [for 3 calls per CNV cluster, kept max cn relative to null 2]
cnv_concor_par$cn %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software),
              FUN = function(vector) {
                  if(min(vector) > 2) {return(max(vector))}
                  else if(max(vector) < 2) {return(min(vector))}
              }) %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "cn" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- cn),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = cn.PennCNV, y = cn.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of copy number status for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(0, 4) +
    ylim(0, 4) +
    theme_minimal()

# Validity assessment - no. of SNPs/probes concordance
cnv_concor_par$probes %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software), FUN = "max") %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "probes" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- probes),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    select(cnv_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "cnv_id",
                names_from = "software",
                values_from = c("cn", "probes", "length", "start_pos", "end_pos"),
                names_sep = ".") %>% 
    ggplot(aes(x = probes.PennCNV, y = probes.QuantiSNP)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5, se = FALSE) +
    geom_abline(slope = 1, color = "red", linetype = "dashed") +
    coord_equal() +
    ggtitle("Concordance of no. of SNPs/probes for partially concordant CNV calls") +
    xlab("PennCNV") +
    ylab("QuantiSNP") +
    xlim(2, 11) +
    ylim(2, 11) +
    theme_minimal()

# Filter harmonized partially concordant calls [Pending: `NA` values filling in `probes`]
Mode <- function(x) {
    ux <- na.omit(unique(x) )
    tab <- tabulate(match(x, ux)); ux[tab == max(tab) ]
}
cnv_concor_par <- cnv_concor_par$cn %>% 
    aggregate(by = list(cnv_concor_par$cnv_id, cnv_concor_par$software),
              FUN = function(vector) {
                  if(min(vector) > 2) {return(max(vector))}
                  else if(max(vector) < 2) {return(min(vector))}
              }) %>% 
    as.data.frame() %>% 
    rename("cnv_id" = "Group.1", "software" = "Group.2", "cn" = "x") %>% 
    left_join(cnv_concor_par %>% 
                  select(- cn),
              by = c("cnv_id", "software")) %>% 
    distinct(cnv_id, software, .keep_all = TRUE) %>% 
    group_by(cnv_id) %>% 
    mutate(inter_length = inter_end_pos - inter_start_pos + 1,
           inter_start_snp = ifelse(inter_start_pos == start_pos, start_snp, as.character(NA)),
           inter_end_snp = ifelse(inter_end_pos == end_pos, end_snp, as.character(NA)),
           inter_probes = ifelse(inter_start_pos == start_pos & inter_end_pos == end_pos, probes, as.numeric(NA))) %>% 
    mutate(across(contains("inter"), ~ Mode(.x))) %>% 
    mutate(across(contains("BF"), ~ Mode(.x))) %>% 
    as.data.frame() %>% 
    select(sample, chr, inter_start_snp, inter_end_snp, inter_probes, cn, inter_length, inter_start_pos, inter_end_pos, contains("BF"), software) %>% 
    pivot_wider(names_from = "software", values_from = "cn") %>% 
    rename_with(~ str_remove(., 'inter_')) %>% 
    rename("cn.penn" = "PennCNV", "cn.quan" = "QuantiSNP")

# Merge completely and partially concordant calls
cnv_concor <- cnv_concor_com %>% mutate(chr = as.numeric(chr)) %>% bind_rows(cnv_concor_par)

# Venn diagram
cnv_venn <- cnv_combine %>% 
    mutate(penn = if_else(software.penn == "PennCNV", 1, 0),
           quan = if_else(software.quan == "QuantiSNP", 1, 0)) %>%
    unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% 
    mutate(penn = if_else(cnv %in% cnv_concor_par_cnvid, 1, penn),
           quan = if_else(cnv %in% cnv_concor_par_cnvid, 1, quan))

cnv_venn <- cnv_venn %>% 
    filter(is.na(penn) | is.na(quan)) %>% 
    select(cnv, penn, quan) %>% 
    bind_rows(cnv_concor %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% mutate(penn = 1, quan = 1))

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))
