## Written by: Kate Johnson

id_list = c("virus",'rep',"AF","allele_freq",
                    "frac","seq_depth","aligner",
                    "tool","parameter",
                    "parameter2")

cov_id_list = c("virus",'rep',"AF","allele_freq",
                    "frac","seq_depth","aligner")

caller_colors = c('#CC2F42', # freebayes - red
                 '#35B6E9', # hc - bright blue
                 '#009E73', # ivar - green
                 '#F9BECD', # lofreq - pink
                 '#F0E442', # mutect2 - yellow
                 '#CCD2EB', # timo - light purple
                 '#9B9E9A', # varscan - grey
                 '#555555') # varscan, custom1 - darker grey)

callers = c('freebayes','hc','ivar',
            'lofreq','mutect2','timo','varscan',
            'varscan_custom1')
names(caller_colors) = callers
caller_colScale_fill <- scale_fill_manual(name = "grp",values = caller_colors)
caller_colScale <- scale_colour_manual(name = "grp",values = caller_colors)



cat_colors  = c('#2F142B', # TP - dark purple
                '#F2A07B', # FP - peach
                '#43668B') # FN - blue
categories = c('TP','FP','FN')
names(cat_colors) = categories
cat_colScale_fill <- scale_fill_manual(name = "grp",values = cat_colors)
cat_colScale <- scale_colour_manual(name = "grp",values = cat_colors)


subsamp_colors = c('#DDDDDD', # 100X - lightest grey
                      '#BBBBBB', # 200X
                      '#999999', # 300X
                      '#777777', # 500X
                      '#555555', # 1000X
                      '#333333', # 10000X
                      '#111111') # 100000X - darkest grey

subsampling_depths = c(100, 200, 300, 500, 1000, 10000, 100000)
names(subsamp_colors) = subsampling_depths
subsamp_colScale_fill <- scale_fill_manual(name = "Coverage",values = subsamp_colors)
subsamp_colScale <- scale_colour_manual(name = "Coverage",values = subsamp_colors)


virus_colors = c('#2F142B', # H1N1 - dark purple
                  '#F2A07B', # H3N2 - peach
                  '#43668B', # VICT - blue
                  '#8DA1E2') # SARS - light purple
strains = c("H1N1","H3N2","VICT","SARS")
names(virus_colors) = strains
strain_colScale_fill <- scale_fill_manual(name = "grp",values = virus_colors)
strain_colScale <- scale_colour_manual(name = "grp",values = virus_colors)

CHROMS = c("PB2","PB1","PA","HA","NP","NA","MP","NS","SARS")

keep_tool_params = c('freebayes_default_NONE',
                    'freebayes_standard_NONE',
                    'freebayes_custom_NONE',
                    'hc_default_NONE',
                    'hc_standard_NONE',
                    'ivar_custom_NONE',
                    'ivar_default_NONE',
                    'ivar_standard_NONE',
                    'lofreq_default_NONE',
                    'lofreq_standard_NONE',
                    'mutect2_custom_unfiltered',
                    'mutect2_default_unfiltered',
                    'mutect2_standard_unfiltered',
                    #'timo_custom_no-binom-check', #######!!!!
                    'timo_custom_NONE', #####!!!!!!!
                    'timo_default_NONE',
                    #'timo_standard_no-binom-check', #######!!!!!!
                    'timo_standard_NONE',  #####!!!!!!!
                    'varscan_custom-1_NONE',
                    'varscan_default_NONE',
                    'varscan_standard_NONE')


all_tool_params = c('freebayes_default_NONE',
                    'freebayes_standard_NONE',
                    'freebayes_custom_NONE',

                    'hc_default_NONE',
                    'hc_standard_NONE',
                    'hc_custom_NONE',

                    'ivar_custom_NONE',
                    'ivar_default_NONE',
                    'ivar_standard_NONE',

                    'lofreq_default_NONE',
                    'lofreq_standard_NONE',
                    'lofreq_custom_NONE',

                    'mutect2_custom_unfiltered',
                    'mutect2_default_unfiltered',
                    'mutect2_standard_unfiltered',

                    'timo_custom_NONE', #####!!!!!!!
                    'timo_default_NONE',
                    'timo_standard_NONE',  #####!!!!!!!
                    'timo_standard_no-binom-check', #######!!!!!!
                    'timo_custom_no-binom-check', #######!!!!

                    'varscan_custom-1_NONE',
                    'varscan_default_NONE',
                    'varscan_standard_NONE')


adj_caller_colors = c('#CC2F42', # freebayes - red
                      '#CC2F42', # freebayes - red
                      '#CC2F42', # freebayes - red

                 '#35B6E9', # hc - bright blue
                 '#35B6E9', # hc - bright blue
                 '#35B6E9', # hc - bright blue

                 '#009E73', # ivar - green
                 '#009E73', # ivar - green
                 '#009E73', # ivar - green

                 '#F9BECD', # lofreq - pink
                 '#F9BECD', # lofreq - pink
                 '#F9BECD', # lofreq - pink

                 '#F0E442', # mutect2 - yellow
                 '#F0E442', # mutect2 - yellow
                 '#F0E442', # mutect2 - yellow


                 '#CCD2EB', # timo - light purple
                 '#CCD2EB', # timo - light purple
                 '#CCD2EB', # timo - light purple

                 '#8DA1E2', # timo-nobino darker purple
                 '#8DA1E2', # timo-nobino darker purple

                 #'#555555', #custom1 - darker gray
                 '#9B9E9A',
                 '#9B9E9A', # varscan - grey
                 '#9B9E9A')


names(adj_caller_colors) = all_tool_params
adj_caller_colScale_fill <- scale_fill_manual(name = "grp",values = adj_caller_colors)
adj_caller_colScale <- scale_colour_manual(name = "grp",values = adj_caller_colors)

################################################################################

rearrange_output_data = function(df, readsim_cov = 100000, id_list, chromlist){

    select_list = c('chrom','pos','af_golden','af_workflow','ref','alt','dp','strain',
                      id_list)

    ntlist = c('A','G','C','T', '-')

    AV = separate(df, sample_id, sep = "_",
               into = c(all_of(id_list))) %>%
                filter(ref %in% ntlist & alt %in% ntlist) %>%
            select(all_of(select_list)) %>%
            droplevels()

    AV$seq_depth = as.numeric(AV$seq_depth)
    AV$coverage = (AV$seq_depth * readsim_cov) #variable called in pipeline, corresponds to highest coverage
    AV$coverage = as.character(AV$coverage)
    AV$CHROM = AV$chrom
    AV = AV %>% separate(CHROM, sep = '_', c("virus_type", 'segment'))
    AV$segment = factor(AV$segment, levels = chromlist)
    AV$parameter2[is.na(AV$parameter2)] = 'NONE'  # have to add so other na's can be turned to zero
    return(AV)
}


add_category_information = function(df){

    return(df %>% mutate(category = ifelse(af_golden == 0 & af_workflow != 0, "FP",
          ifelse(af_golden != 0 & af_workflow == 0, "FN",
                 ifelse(af_golden != 0 & af_workflow != 0, "TP", "OTHER")))))
}


rearrange_coverage = function(df, cov_id_list){
  df = as.data.frame(df)
  df$sample_id = df$name
  df = separate(df, sample_id, sep = "_",
             into = c(all_of(cov_id_list))) %>%
        unique() %>%
        droplevels()
  df$seq_depth = as.numeric(df$seq_depth)
  return(df)
}

adjust_counts = function(count_df){

    count_df = count_df %>%
                pivot_wider(names_from = category, values_from = category_count)

    if (!"FP" %in% colnames(count_df)){
      print("No false-positives - adjusting")
      count_df$FP = 0

    }else{print("FP present - no adjustment needed")}
    count_df[is.na(count_df)] = 0
    count_df$TPR = count_df$TP/(count_df$TP + count_df$FN)
    count_df$PPV = count_df$TP/(count_df$TP + count_df$FP)
    count_df$FNR = count_df$FN/(count_df$FN + count_df$TP)
    count_df$FDR = count_df$FP/(count_df$FP + count_df$TP)
    count_df$F1 = 2 * ((count_df$PPV * count_df$TPR)/(count_df$PPV + count_df$TPR))
    count_df[is.na(count_df)] = 0
    count_df$coverage = as.numeric(as.character(count_df$coverage))
    return(count_df)
}


adjust_information_types = function(df){
    df$coverage = as.numeric(as.character(df$coverage))
    df$tool_param = paste0(df$tool, '_', df$parameter)
    df$total_input_snps = df$TP + df$FN
    df$tool_strain = paste0(df$strain, '_', df$tool)
    df$total_positives = df$FP + df$TP
    df$total_snps = df$FP + df$TP + df$FN
    df$tp_ratio = df$TP/df$total_snps
    df$fp_ratio = df$FP/df$total_snps
    df$fn_ratio = df$FN/df$total_snps
    return(df)
}


PlotTheme1 = theme_bw() +
              theme(legend.key = element_blank(),
                strip.background = element_rect(colour="black", fill="white"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

PlotTheme2 = theme_bw() +
                theme(legend.key = element_blank(),
                strip.background = element_rect(colour="black", fill="white"))

PlotTheme3 = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.key = element_blank(),
                strip.background = element_rect(colour="black", fill="white"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

PlotTheme4 = theme(  panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.key = element_blank(),
                      strip.background = element_rect(colour="black", fill="white"),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
