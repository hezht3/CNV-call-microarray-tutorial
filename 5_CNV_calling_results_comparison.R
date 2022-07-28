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
# Rules:
# Group by sample, chromosome: sort CNV calls by starting position (first), ending position (second)
# Calculate difference of starting position of current call with ending position of previous call AND
#           difference of ending position of current call with starting position of previous call
# For partially concordant calls, product of the above 2 difference < 0

# Illustrations:
# Rest CNV calls can be partitioned in to following scenairos:
# + ----------------------------------------------------------------------------------------------------------------------- +
# + - CNV call visulization - | - Start pos - lag(end pos) - | - End pos - lag(start pos) - | - Product of left 2 columns - |
# + ----------------------------------------------------------------------------------------------------------------------- +
# + ----------------------------------------------------------------------------------------------------------------------- +
# | Partially concordant                                                                                                    |
# + ----------------------------------------------------------------------------------------------------------------------- +
# | ======================    |              < 0             |             > 0              |             < 0               |
# | =========                 |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# |              =========    |              < 0             |             > 0              |             < 0               |
# | ======================    |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# | ======================    |              < 0             |             > 0              |             < 0               |
# |       ==========          |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# |       ==========          |              < 0             |             > 0              |             < 0               |
# | ======================    |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# | ===============           |              < 0             |             > 0              |             < 0               |
# |        ===============    |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# |        ===============    |              < 0             |             > 0              |             < 0               |
# | ===============           |                              |                              |                               |
# + ----------------------------------------------------------------------------------------------------------------------- +
# + Not concordant                                                                                                          |
# + ----------------------------------------------------------------------------------------------------------------------- +
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# |                  =======  |              < 0             |             < 0              |             > 0               |
# | ========                  |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +
# | ========                  |              > 0             |             > 0              |             > 0               |
# |                  =======  |                              |                              |                               |
# + ------------------------- + ---------------------------- + ---------------------------- + ----------------------------- +

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
    mutate(id_num = 1:n()) %>% 
    mutate(start_end_diff = start_pos - lag(end_pos),
           end_start_diff = end_pos - lag(start_pos)) %>% 
    mutate(par_indi = if_else(start_end_diff * end_start_diff < 0, id_num, as.integer(NA))) %>% 
    mutate(par_indi = if_else(!is.na(par_indi) & !is.na(lag(par_indi)), lag(par_indi), par_indi)) %>% 
    mutate(par_indi = if_else(!is.na(lead(start_end_diff)) & !is.na(lead(end_start_diff)) & lead(start_end_diff) * lead(end_start_diff) < 0,
                              lead(par_indi), par_indi)) %>% 
    filter(!is.na(par_indi)) %>% 
    mutate(par_pair_id = paste(sample, chr, par_indi, sep = "_"))
cnv_concor_par_cnvid <- cnv_concor_par %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% pull(cnv)

# Intersection [Pending resolve: finding intersection for partially concordant calls >= 3]
cnv_concor_par <- cnv_concor_par %>% 
    group_by(par_pair_id) %>% 
    mutate(id_num = 1:n()) %>% 
    mutate(inter_start_pos = max(start_pos), inter_end_pos = min(end_pos)) %>% 
    mutate(inter_start_pos = case_when(n() == 3 & id_num == 1 & start_pos >= lead(start_pos) ~ start_pos,
                                       n() == 3 & id_num == 1 & start_pos < lead(start_pos) ~ lead(start_pos),
                                       n() == 3 & id_num == 3 & start_pos >= lag(start_pos) ~ start_pos,
                                       n() == 3 & id_num == 3 & start_pos < lag(start_pos) ~ lag(start_pos),
                                       TRUE ~ inter_start_pos),
           inter_end_pos = case_when(n() == 3 & id_num == 1 & end_pos >= lead(end_pos) ~ lead(end_pos),
                                     n() == 3 & id_num == 1 & end_pos < lead(end_pos) ~ end_pos,
                                     n() == 3 & id_num == 3 & end_pos >= lag(end_pos) ~ end_pos,
                                     n() == 3 & id_num == 3 & end_pos < lag(end_pos) ~ lag(end_pos),
                                     TRUE ~ inter_end_pos)) %>% 
    mutate(inter_length = inter_end_pos - inter_start_pos + 1) %>% 
    mutate(inter_length = if_else(n() == 3 & id_num == 2, lag(inter_length) + lead(inter_length), inter_length))

# Evaluate intersection regions - visualization [Pending resolve: finding intersection for partially concordant calls >= 3]
cnv_concor_par %>% distinct(par_pair_id) %>% pull(par_pair_id) %>% 
    map(~ cnv_concor_par %>% 
            filter(par_pair_id == .x) %>% 
            mutate(cnv_call = paste(sample, chr, "call", id_num, software, sep = ".")) %>% 
            ggplot(aes(y = cnv_call, color = software)) +
            geom_linerange(aes(xmin = start_pos, xmax = end_pos), size = 3) +
            xlab("Position") +
            ylab("Sample ID.CHR.CNV.Software") +
            scale_color_manual("Software", values = c("#ff9900", "#146eb4")) +
            theme_minimal())

# Evaluate intersection regions - percent agreement [Pending resolve: finding intersection for partially concordant calls >= 3]
cnv_concor_par %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
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

cnv_concor_par %>% 
    mutate(perc.agree = if_else(inter_length / length <= 1, inter_length / length, 1)) %>% 
    mutate(software = factor(software,
                             levels = c("PennCNV", "QuantiSNP"),
                             labels = c("Gold standard: PennCNV", "Gold standard: QuantiSNP"))) %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    ggplot() +
    geom_histogram(aes(x = perc.agree, fill = software)) +
    facet_wrap(. ~ software, scales = "free_x") +
    geom_vline(aes(xintercept = 0.6, color = "red"), linetype = "dashed", show.legend = F) +
    xlab("Percent agreement") +
    ylab("Frequency") +
    scale_fill_manual("Software", values = c("#ff9900", "#146eb4"), guide = "none") +
    theme_minimal()

# Validity assessment - length concordance
cnv_concor_par %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    select(par_pair_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "par_pair_id",
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

# Validity assessment - start position concordance
cnv_concor_par %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    select(par_pair_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "par_pair_id",
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

# Validity assessment - end position concordance
cnv_concor_par %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    select(par_pair_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "par_pair_id",
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

# Validity assessment - copy number
cnv_concor_par %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    select(par_pair_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "par_pair_id",
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
cnv_concor_par %>% 
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    select(par_pair_id, cn, probes, length, software, start_pos, end_pos) %>% 
    pivot_wider(id_cols = "par_pair_id",
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

# Filter harmonized partially concordant calls
Mode <- function(x) {
    ux <- na.omit(unique(x) )
    tab <- tabulate(match(x, ux)); ux[tab == max(tab) ]
}
cnv_concor_par <- cnv_concor_par %>%
    filter(!(sample == "99HI0698C" & chr == 14 & probes == 9)) %>% 
    group_by(par_pair_id) %>% 
    mutate(inter_start_snp = ifelse(inter_start_pos == start_pos, start_snp, as.character(NA)),
           inter_end_snp = ifelse(inter_end_pos == end_pos, end_snp, as.character(NA)),
           inter_probes = ifelse(inter_start_pos == start_pos & inter_end_pos == end_pos, probes, as.numeric(NA))) %>% 
    mutate(inter_start_pos = ifelse(par_pair_id == "99HI0698C_14_4", start_pos[software == "PennCNV"], inter_start_pos),
           inter_end_pos = ifelse(par_pair_id == "99HI0698C_14_4", end_pos[software == "PennCNV"], inter_end_pos),
           inter_start_snp = ifelse(par_pair_id == "99HI0698C_14_4", start_snp[software == "PennCNV"], inter_start_snp),
           inter_end_snp = ifelse(par_pair_id == "99HI0698C_14_4", end_snp[software == "PennCNV"], inter_end_snp),
           inter_probes = ifelse(par_pair_id == "99HI0698C_14_4", probes[software == "PennCNV"], inter_probes),
           inter_length = ifelse(par_pair_id == "99HI0698C_14_4",
                                 end_pos[software == "PennCNV"] - start_pos[software == "PennCNV"] + 1,
                                 inter_length)) %>% 
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
    filter(!(cnv %in% cnv_concor_par_cnvid)) %>% 
    bind_rows(cnv_venn %>% filter(cnv %in% (cnv_concor_par %>% unite("cnv", c(sample, chr, start_snp, end_snp), sep = ".") %>% pull(cnv))))

ggvenn(list(`PennCNV` = cnv_venn %>% filter(penn == 1) %>% pull(cnv),
            `QuantiSNP` = cnv_venn %>% filter(quan == 1) %>% pull(cnv)))
