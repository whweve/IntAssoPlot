#' @title plot association with LD and annotation at a given gene
#' @description  This function plot the association with
#' linkage disequiblism and annotation at the level of a single gene.
#' @author Hongwei Wang <\email{whweve@163.com}>
#' @param transcript the transcript of gene, required.
#' @param gtf the annotation file, required.
#' @param association the association table, required.
#' @param hapmap the genotype file for computing leadsnpLD in the format of hapmap. The file should be the same file used for coumputing association. required. 
#' @param hapmap_ld the genotype file for computing trangleLD in the format of hapmap, not required. If hapmap_ld was not provided, hapmap would be used.
#' @param slide_length the sliding window length for computing LD, default -1.
#' @param threadN the number of (CPU) cores used for computing LD, default 1.
#' @param up the upper distance from the start position of gene
#' @param down the down distance from the end position of gene
#' @param threshold the significant level of the assocition, default NULL.
#' @param ldstatistics the statistics used for computing LD, default rsquare, and the optional is dprime.
#' @param leadsnp snp name provided by user
#' @param link2gene a dataframe speicify markers to be linked from GWAS to genic structure, default NULL. When link2gene is 'NULL', locis that passed the threshold will be linked. Please see help('marker2link').
#' @param triangleLD show LD in the format lile triangle, default TRUE.
#' @param link2LD a dataframe speicify markers to be linked from genic structure to LD matrix, default NULL. When link2gene is 'NULL', locis that passed the threshold will be linked. Please see help('marker2link').
#' @param leadsnpLD show LD of the locis when compared with the most significant loci, default TRUE.
#' @param exon_colour the colour of exon, default gray.
#' @param cds_colour the colour of cds, default black.
#' @param utr_colour the colour of utr, default gray.
#' @param intron_colour the colour of intron, default gray.
#' @param colour02 the colour of LD statistics ranged between 0.0 and 0.2, default gray.
#' @param colour04 the colour of LD statistics ranged between 0.2 and 0.4, default cyan.
#' @param colour06 the colour of LD statistics ranged between 0.4 and 0.6, default green.
#' @param colour08 the colour of LD statistics ranged between 0.6 and 0.8, default yellow.
#' @param colour10 the colour of LD statistics ranged between 0.8 and 1.0, default red.
#' @param leadsnp_shape the shape of leadsnp, default 23. For others, please see help('points').
#' @param leadsnp_colour the shape of the point of leadsnp, default black. For others, please see help('points').
#' @param leadsnp_fill the filled colour of the point of leadsnp, default purple. For others, please see help('points').
#' @param leadsnp_size the size of of the point of leadsnp, default 1.5. For others, please see help('points').
#' @param marker2highlight a dataframe speicify markers to be showed by colour,shape,fill,size, default NULL. Please see help('marker2highlight').
#' @param marker2label a dataframe speicify markers to be labeled, default NULL. Please see help('marker2link')
#' @param marker2label_angle angel of labeled text, default 60.
#' @param marker2label_size size of labeled text, default 1.
#' @param thresholdlinecolour colour of threshold line, default gray.
#' @param upperpointsize size of point of association sites, default 1.
#' @return ggplot2 plot
#' @export
#' @import ggplot2 SNPRelate reshape2 gdsfmt ggrepel
#' @examples 
#' data(gtf)
#' data(zmvpp1_association)
#' data(zmvpp1_hapmap)
#' data(marker2highlight)
#' data(marker2link)
#' IntGenicPlot('GRMZM2G170927_T01',gtf,association=zmvpp1_association,hapmap=zmvpp1_hapmap,
#' hapmap_ld = zmvpp1_hapmap,threshold=8,up=500,down=600,leadsnpLD = FALSE,
#' marker2highlight=marker2highlight,link2gene=marker2link,link2LD=marker2link,
#' marker2label=marker2link,marker2label_angle=60,marker2label_size=2)
IntGenicPlot <- function(transcript, gtf, association, hapmap, hapmap_ld = NULL, 
                         slide_length = -1, threadN = 1, up = NULL, down = NULL, threshold = NULL, ldstatistics = "rsquare", 
                         leadsnp = NULL, link2gene = NULL, triangleLD = TRUE, link2LD = NULL, leadsnpLD = TRUE, 
                         exon_colour = "gray", cds_colour = "black", utr_colour = "gray", intron_colour = "gray", 
                         colour02 = "gray", colour04 = "cyan", colour06 = "green", colour08 = "yellow", 
                         colour10 = "red", leadsnp_shape = 23, leadsnp_colour = "black", leadsnp_fill = "purple", 
                         leadsnp_size = 1.5, marker2highlight = NULL, marker2label = NULL, marker2label_angle = 60, 
                         marker2label_size = 1,thresholdlinecolour="gray",upperpointsize=1) {
  if (sum(grepl(transcript, gtf$V9)) == 0) {
    stop("please provide the correct transcript or the gtf file")
  } else {
    
    if (names(association) %in% c("Marker", "Locus",  "Site",   "p") %>% sum() != 4 ) {
      colpos <- which(!(c("Marker", "Locus",  "Site",   "p") %in% names(association)))
      print(paste0("Required column ",c("Marker", "Locus",  "Site",   "p")[colpos]," not existed, this may lead to error. See IntAssoPlot::association for help"))
    } 
    if (names(association) %in% c("Marker", "Locus",  "Site",   "p") %>% sum() == 4 ) {
      print("checking association table. Done")
    }
    
    if (grepl("gene_id (\\S+) .+",gtf$V9) %>% sum() >=1) {
      print("checking gtf table. Done")
    } 
    if (grepl("gene_id (\\S+) .+",gtf$V9) %>% sum() == 0) {
      print("No text like 'gene_id xxxxx' appeared in the ninth column of gtf. This may lead to error. Rebuild the gtf file from gff file using gffread. See IntAssoPlot::gtf for help")
    }
    
    if(names(hapmap)[1] != "rs") {
    print("Converting the first column name to rs")
    names(hapmap)[1] = "rs"
    }
    if(sum(asso$Marker %in% hapmap$rs) >= 1) {
      print("There are identical marker names within association file and the hapmap file. Done")
    }
    
    if( asso$Marker %in% hapmap$rs %>% sum() ==0 ) {
      print("There are no identical marker names within association file and the hapmap file. This may lead to error")
    }
    
    R2 <- Site <- Site2 <- V4 <- V5 <- V9 <- group <- p <- NULL
    transcript_structure_cds_list <- transcript_structure_exon_list <- NULL
    transcript_structure_utr_list <- x <- xend <- y <- yend <- NULL
    # globalVariables(names(gtf)) globalVariables(names(association))
    # globalVariables(names(hapmap))
    transcript_corrdination <- gtf[grepl(transcript, gtf$V9), ]
    chromosome_association <- association[association$Locus == unique(transcript_corrdination$V1), 
    ]
    transcript_corrdination <- gtf[grepl(transcript, gtf$V9), ]
    transcript_min <- ifelse(is.null(up), min(transcript_corrdination$V4), 
                             min(transcript_corrdination$V4) - up)
    transcript_max <- ifelse(is.null(down), max(transcript_corrdination$V5), 
                             max(transcript_corrdination$V5) + down)
    transcript_association <- chromosome_association[chromosome_association$Site >= 
                                                       transcript_min & chromosome_association$Site <= transcript_max, ]
    transcript_association <- transcript_association[order(transcript_association$Site), 
    ]
    if (dim(transcript_association)[1] < 2) {
      stop("Less than 2 markers, can not compute LD.")
    } else {
      # compute the meta variable
      pvalue_range <- pretty(c(0,-log10(transcript_association$p)))
      # adjust the yaxis to fit in the LD plot
      fold <- ((transcript_max - transcript_min) * 2/3)/max(pvalue_range)
      n_pvalue_range <- length(pvalue_range)
      marker_number = dim(transcript_association)[1]
      length = (transcript_max - transcript_min)
      distance = 0.5 * length/(marker_number - 1)
      
      # transcript start and end
      for (struct in c("utr", "cds", "exon")) {
        assign(paste0("transcript_structure_", struct), transcript_corrdination[grep(struct, 
                                                                                     transcript_corrdination$V3, ignore.case = TRUE), ])
        if (dim(get(paste0("transcript_structure_", struct)))[1] > 0) {
          assign(paste0("transcript_structure_", struct, "_list"), list(geom_segment(data = get(paste0("transcript_structure_", 
                                                                                                       struct)), aes(x = V4, xend = V5, y = -max(pvalue_range) * fold/30, 
                                                                                                                     yend = -max(pvalue_range) * fold/30), colour = get(paste0(struct, 
                                                                                                                                                                               "_colour")), size = 4)))
        } else {
          assign(paste0("transcript_structure_", struct, "_list"), NULL)
        }
      }
      transcript_intron_structure <- list(geom_segment(aes(x = transcript_min, 
                                                           xend = transcript_max, y = -max(pvalue_range) * fold/30, yend = -max(pvalue_range) * 
                                                             fold/30), colour = intron_colour, size = 1))
      # decide whether to rotate x axis
      scale_x <- ifelse(unique(transcript_corrdination$V7) == "-", list(scale_x_reverse(limits = c((transcript_max - 
                                                                                                      transcript_min)/6 + transcript_max, transcript_min), breaks = seq(transcript_max, 
                                                                                                                                                                        transcript_min, transcript_min - transcript_max))), list(scale_x_continuous(limits = c(transcript_min - 
                                                                                                                                                                                                                                                                 (transcript_max - transcript_min)/6, transcript_max), breaks = seq(transcript_min, 
                                                                                                                                                                                                                                                                                                                                    transcript_max, transcript_max - transcript_min))))
      # label the yaxis
      scale_y_line <- ifelse(unique(transcript_corrdination$V7) == "-", list(geom_segment(aes(x = (transcript_max - 
                                                                                                     transcript_min)/30 + transcript_max, y = min(pvalue_range), xend = (transcript_max - 
                                                                                                                                                                           transcript_min)/30 + transcript_max, yend = max(pvalue_range) * 
                                                                                                fold))), list(geom_segment(aes(x = transcript_min - (transcript_max - 
                                                                                                                                                       transcript_min)/30, y = min(pvalue_range), xend = transcript_min - 
                                                                                                                                 (transcript_max - transcript_min)/30, yend = max(pvalue_range) * 
                                                                                                                                 fold))))
      scale_y_ticks <- ifelse(unique(transcript_corrdination$V7) == "-", 
                              list(geom_segment(aes(x = rep((transcript_max - transcript_min)/15 + 
                                                              transcript_max, n_pvalue_range), y = pvalue_range * fold, xend = (transcript_max - 
                                                                                                                                  transcript_min)/30 + transcript_max, yend = pvalue_range * fold))), 
                              list(geom_segment(aes(x = rep(transcript_min - (transcript_max - 
                                                                                transcript_min)/15, n_pvalue_range), y = pvalue_range * fold, 
                                                    xend = rep(transcript_min - (transcript_max - transcript_min)/30, 
                                                               n_pvalue_range), yend = pvalue_range * fold))))
      scale_y_text <- ifelse(unique(transcript_corrdination$V7) == "-", list(geom_text(aes(x = rep((transcript_max - 
                                                                                                      transcript_min)/12 + transcript_max, n_pvalue_range), y = pvalue_range * 
                                                                                             fold, label = pvalue_range))), list(geom_text(aes(x = rep(transcript_min - 
                                                                                                                                                         (transcript_max - transcript_min)/12, n_pvalue_range), y = pvalue_range * 
                                                                                                                                                 fold, label = pvalue_range))))
      # add threshold line
      if (is.null(threshold)) {
        threshold_line <- list(NULL)
      }
      if (all(length(threshold) > 0, threshold <= max(pvalue_range))) {
        threshold_line <- list(geom_segment(aes(x = transcript_min, xend = transcript_max, 
                                                y = threshold * fold, yend = threshold * fold), linetype = "longdash", 
                                            colour = thresholdlinecolour))
      }
      if (all(length(threshold) > 0, threshold > max(pvalue_range))) {
        threshold_line <- list(NULL)
        print("no -log10(p) pass the threshold, will not draw threshold line")
      }
      # compute the LD, leadsnp or triangle
      if (any(isTRUE(leadsnpLD), isTRUE(triangleLD)) & is.null(hapmap)) {
        print("no hapmap data found, please provide the hapmap")
        ld_leadsnp_colour <- list(NULL)
        bottom_trianglLD = list(NULL)
      }
      if (all(!isTRUE(leadsnpLD), !isTRUE(triangleLD), !is.null(hapmap))) {
        ld_leadsnp_colour <- list(NULL)
        bottom_trianglLD <- list(NULL)
      }
      # link association and LD for the significant loci link between LD and genic
      # structure
      if (any(isTRUE(leadsnpLD), isTRUE(triangleLD)) & !is.null(hapmap)) {
        names(hapmap) <- sub("#", "", names(hapmap))
        gene_snp <- hapmap[hapmap$chrom == unique(transcript_corrdination$V1) & 
                             hapmap$pos >= transcript_min & hapmap$pos <= transcript_max, 
        ]
        # gene_snp = hapmap2
        names(gene_snp) <- sub("#", "", names(gene_snp))
        gene_snp <- gene_snp[!duplicated(gene_snp$rs), ]
        # convert the SNP to numeric format
        major_allele <- paste0(substr(gene_snp$allele, 1, 1), substr(gene_snp$allele, 
                                                                     1, 1))
        minor_allele <- paste0(substr(gene_snp$allele, 3, 3), substr(gene_snp$allele, 
                                                                     3, 3))
        heter_left <- paste0(substr(gene_snp$allele, 1, 1), substr(gene_snp$allele, 
                                                                   3, 3))
        heter_right <- paste0(substr(gene_snp$allele, 3, 3), substr(gene_snp$allele, 
                                                                    1, 1))
        # if allele equal to major allele, 0, else 2
        for (j in 12:dim(gene_snp)[2]) {
          gene_snp[gene_snp[, j] == major_allele, j] = 2
          gene_snp[gene_snp[, j] == minor_allele, j] = 0
          gene_snp[gene_snp[, j] == "NN", j] = NA
          heter_position_left <- which(isTRUE(gene_snp[, j] == heter_left))
          heter_position_right <- which(isTRUE(gene_snp[, j] == heter_right))
          if (length(heter_position_left) > 1) {
            gene_snp[heter_position_left, j] = 1
          }
          if (length(heter_position_right) > 1) {
            gene_snp[heter_position_right, j] = 1
          }
        }
        gene_snp2 <- gene_snp[, 12:dim(gene_snp)[2]]
        gene_snp2 <- as.matrix(sapply(gene_snp2, as.numeric))
        
        snpgdsCreateGeno("test.gds", genmat = gene_snp2, sample.id = names(gene_snp)[12:dim(gene_snp)[2]], 
                         snp.id = gene_snp$rs, snp.position = gene_snp$pos, snp.allele = gene_snp$alleles, 
                         snpfirstdim = TRUE)
        
        genofile <- snpgdsOpen("test.gds")
        if (ldstatistics == "rsquare") {
          aa = snpgdsLDMat(genofile, slide = slide_length, method = "corr", 
                           num.thread = threadN)
        }
        if (ldstatistics == "dprime") {
          aa = snpgdsLDMat(genofile, slide = slide_length, method = "dprime", 
                           num.thread = threadN)
        }
        snpgdsClose(genofile)
        ld = aa$LD
        if (ldstatistics == "rsquare") 
          ld <- ld^2
        names(ld) <- gene_snp$rs
        ld <- melt(ld)
        
        marker_info <- data.frame(index = 1:dim(gene_snp)[1], marker_name = gene_snp$rs)
        ld$Var1 <- marker_info$marker[match(ld$Var1, marker_info$index)]
        ld$Var2 <- marker_info$marker[match(ld$Var2, marker_info$index)]
        if (ldstatistics == "rsquare") {
          lengend_name = expression(italic(r)^2)
        } else if (ldstatistics == "dprime") {
          lengend_name = expression(D * {
            "'"
          })
        }
        ld <- ld[!is.na(ld$value), ]
        ld_reverse <- data.frame(Var1 = ld$Var2, Var2 = ld$Var1, value = ld$value)
        ld <- rbind(ld, ld_reverse)
        marker_pos <- transcript_association[, c("Marker", "Site")]
        ld$Site1 <- marker_pos$Site[match(ld$Var1, marker_pos$Marker)]
        ld$Site2 <- marker_pos$Site[match(ld$Var2, marker_pos$Marker)]
        # ld <- merge(ld,marker_pos,by.x='Var1',by.y = 'Marker') ld <-
        # merge(ld,marker_pos,by.x='Var2',by.y = 'Marker') names(ld) =
        # sub('Site.x','Site1',names(ld)) names(ld) = sub('Site.y','Site2',names(ld))
        if (isTRUE(leadsnpLD)) {
          if (is.null(leadsnp)) {
            leadsnp <- as.character(transcript_association[which.min(transcript_association$p), 
                                                           "Marker"])
          }
          if (!is.null(leadsnp)) {
            leadsnp <- leadsnp
          }
          ld_leadsnp <- ld[ld$Var1 == leadsnp, ]
          ld_leadsnp <- merge(ld_leadsnp, transcript_association, by.x = "Var2", 
                              by.y = "Marker")
          ld_leadsnp$R2 <- 0.2 * (ld_leadsnp$value%/%0.2 + as.logical(ld_leadsnp$value%/%0.2))
          ld_leadsnp$R2 <- as.character(ld_leadsnp$R2)
          ld_leadsnp$R2[ld_leadsnp$R2 == "0"] = "0.2"
          ld_leadsnp$R2[ld_leadsnp$R2 == "1.2"] = "1"
          if (length(which(ld_leadsnp$Var1 == leadsnp & ld_leadsnp$Var2 == 
                           leadsnp)) >= 1) {
            ld_leadsnp <- ld_leadsnp[!(ld_leadsnp$Var1 == leadsnp & ld_leadsnp$Var2 == 
                                         leadsnp), ]
          }
          ld_leadsnp_colour <- list(geom_point(data = ld_leadsnp, aes(Site2, 
                                                                      -log10(p) * fold, fill = R2), shape = 21, colour = "black",size=upperpointsize), 
                                    scale_fill_manual(values = c(`0.2` = colour02, `0.4` = colour04, 
                                                                 `0.6` = colour06, `0.8` = colour08, `1` = colour10), labels = c("0-0.2", 
                                                                                                                                 "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), name = lengend_name))
        }
        if (!isTRUE(leadsnpLD)) {
          ld_leadsnp_colour <- list(NULL)
        }
        if (isTRUE(triangleLD)) {
          if (is.null(hapmap_ld)) {
            hapmap_ld <- hapmap
          }
          hapmap_ld <- hapmap_ld[hapmap_ld$chrom == unique(transcript_corrdination$V1) & 
                                   hapmap_ld$pos >= transcript_min & hapmap_ld$pos <= transcript_max, 
          ]
          marker_number = dim(hapmap_ld)[1]
          length = (transcript_max - transcript_min)
          distance = 0.5 * length/(marker_number - 1)
          names(hapmap_ld) <- sub("#", "", names(hapmap_ld))
          gene_snp <- hapmap_ld
          names(gene_snp) <- sub("#", "", names(gene_snp))
          gene_snp <- gene_snp[!duplicated(gene_snp$rs), ]
          # convert the SNP to numeric format
          major_allele <- paste0(substr(gene_snp$allele, 1, 1), substr(gene_snp$allele, 
                                                                       1, 1))
          minor_allele <- paste0(substr(gene_snp$allele, 3, 3), substr(gene_snp$allele, 
                                                                       3, 3))
          heter_left <- paste0(substr(gene_snp$allele, 1, 1), substr(gene_snp$allele, 
                                                                     3, 3))
          heter_right <- paste0(substr(gene_snp$allele, 3, 3), substr(gene_snp$allele, 
                                                                      1, 1))
          # if allele equal to major allele, 0, else 2
          for (j in 12:dim(gene_snp)[2]) {
            gene_snp[gene_snp[, j] == major_allele, j] = 2
            gene_snp[gene_snp[, j] == minor_allele, j] = 0
            gene_snp[gene_snp[, j] == "NN", j] = NA
            heter_position_left <- which(isTRUE(gene_snp[, j] == heter_left))
            heter_position_right <- which(isTRUE(gene_snp[, j] == heter_right))
            if (length(heter_position_left) > 1) {
              gene_snp[heter_position_left, j] = 1
            }
            if (length(heter_position_right) > 1) {
              gene_snp[heter_position_right, j] = 1
            }
          }
          gene_snp2 <- gene_snp[, 12:dim(gene_snp)[2]]
          gene_snp2 <- as.matrix(sapply(gene_snp2, as.numeric))
          snpgdsCreateGeno("test.gds", genmat = gene_snp2, sample.id = names(gene_snp)[12:dim(gene_snp)[2]], 
                           snp.id = gene_snp$rs, snp.position = gene_snp$pos, snp.allele = gene_snp$alleles, 
                           snpfirstdim = TRUE)
          genofile <- snpgdsOpen("test.gds")
          if (ldstatistics == "rsquare") {
            aa = snpgdsLDMat(genofile, slide = slide_length, method = "corr", 
                             num.thread = threadN)
          }
          if (ldstatistics == "dprime") {
            aa = snpgdsLDMat(genofile, slide = slide_length, method = "dprime", 
                             num.thread = threadN)
          }
          snpgdsClose(genofile)
          ld = aa$LD
          if (ldstatistics == "rsquare") 
            ld <- ld^2
          names(ld) <- gene_snp$rs
          ld <- melt(ld)
          marker_info <- data.frame(index = 1:dim(gene_snp)[1], marker_name = gene_snp$rs)
          ld$Var1 <- marker_info$marker[match(ld$Var1, marker_info$index)]
          ld$Var2 <- marker_info$marker[match(ld$Var2, marker_info$index)]
          if (ldstatistics == "rsquare") {
            lengend_name = expression(italic(r)^2)
          } else if (ldstatistics == "dprime") {
            lengend_name = expression(D * {
              "'"
            })
          }
          ld <- ld[!is.na(ld$value), ]
          ld_reverse <- data.frame(Var1 = ld$Var2, Var2 = ld$Var1, value = ld$value)
          ld <- rbind(ld, ld_reverse)
          marker_pos <- hapmap_ld[, c("rs", "pos")]
          ld$Site1 <- marker_pos$pos[match(ld$Var1, marker_pos$rs)]
          ld$Site2 <- marker_pos$pos[match(ld$Var2, marker_pos$rs)]
          # ld <- merge(ld,marker_pos,by.x='Var1',by.y = 'rs') ld <-
          # merge(ld,marker_pos,by.x='Var2',by.y = 'rs') names(ld) =
          # sub('pos.x','Site1',names(ld)) names(ld) = sub('pos.y','Site2',names(ld))
          # compute the LD position, the sequence ranged from small to big
          marker_pair = NULL
          center_x = NULL
          center_y = NULL
          locib <- rep(1:(marker_number - 1), (marker_number - 1):1)
          locia <- NULL
          for (row in 1:(marker_number - 1)) {
            locia <- c(locia, seq(1:(marker_number - row)))
          }
          marker_pair <- paste0(locia, "_", locia + locib)
          center_x <- distance * (locia + locia + locib - 2)
          center_y <- -locib * distance
          upper_center_x <- center_x
          upper_center_y <- center_y + distance
          lower_center_x <- center_x
          lower_center_y <- center_y - distance
          left_center_x <- center_x - distance
          left_center_y <- center_y
          right_center_x <- center_x + distance
          right_center_y <- center_y
          poly_data <- data.frame(group = rep(marker_pair, 4), x = c(upper_center_x, 
                                                                     right_center_x, lower_center_x, left_center_x) + transcript_min, 
                                  y = c(upper_center_y, right_center_y, lower_center_y, left_center_y) - 
                                    4 * max(pvalue_range) * fold/30, label = rep(c(1, 2, 3, 4), 
                                                                                 each = length(upper_center_x)))
          poly_data$marker1 <- sub("([0-9]+)_[0-9]+", "\\1", poly_data$group)
          poly_data$marker2 <- sub("[0-9]+_([0-9]+)", "\\1", poly_data$group)
          # transcript_association <-
          # transcript_association[order(transcript_association$Site),]
          # transcript_association$marker_number <- 1:dim(transcript_association)[1]
          # marker_index <- transcript_association[,c('Marker','marker_number')]
          hapmap_ld <- hapmap_ld[order(hapmap_ld$pos), ]
          hapmap_ld$marker_number <- 1:dim(hapmap_ld)[1]
          marker_index <- hapmap_ld[, c("rs", "marker_number")]
          poly_data$Var1 <- marker_index$rs[match(poly_data$marker1, marker_index$marker_number)]
          poly_data$Var2 <- marker_index$rs[match(poly_data$marker2, marker_index$marker_number)]
          # poly_data <- merge(poly_data,marker_index,by.x='marker1',by.y =
          # 'marker_number') poly_data <-
          # merge(poly_data,marker_index,by.x='marker2',by.y = 'marker_number')
          # names(poly_data) = sub('rs.x','Var1',names(poly_data)) names(poly_data) =
          # sub('rs.y','Var2',names(poly_data))
          poly_data$value <- ld$value[match(paste0(poly_data$Var1, "/", 
                                                   poly_data$Var2), paste0(ld$Var1, "/", ld$Var2))]
          poly_data$Site1 <- ld$Site1[match(paste0(poly_data$Var1, "/", 
                                                   poly_data$Var2), paste0(ld$Var1, "/", ld$Var2))]
          poly_data$Site2 <- ld$Site2[match(paste0(poly_data$Var1, "/", 
                                                   poly_data$Var2), paste0(ld$Var1, "/", ld$Var2))]
          # poly_data <- merge(poly_data,ld,by.x=c('Var1','Var2'),by.y =
          # c('Var1','Var2'))
          poly_data$R2 <- 0.2 * (poly_data$value%/%0.2 + as.logical(poly_data$value%/%0.2))
          poly_data$R2 <- as.character(poly_data$R2)
          poly_data$R2[poly_data$R2 == "0"] = "0.2"
          poly_data$R2[poly_data$R2 == "1.2"] = "1"
          poly_data <- poly_data[order(poly_data$group, poly_data$label), 
          ]
          if (!isTRUE(leadsnpLD)) {
            bottom_trianglLD = list(geom_polygon(data = poly_data, aes(group = group, 
                                                                       x = x, y = y - (transcript_max - transcript_min)/50, fill = R2)), 
                                    scale_fill_manual(values = c(`0.2` = colour02, `0.4` = colour04, 
                                                                 `0.6` = colour06, `0.8` = colour08, `1` = colour10), labels = c("0-0.2", 
                                                                                                                                 "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), name = lengend_name))
          }
          if (isTRUE(leadsnpLD)) {
            bottom_trianglLD = list(geom_polygon(data = poly_data, aes(group = group, 
                                                                       x = x, y = y - (transcript_max - transcript_min)/50, fill = R2)))
          }
          
        }
        if (!isTRUE(triangleLD)) {
          bottom_trianglLD <- list(NULL)
        }
      }
      # link line from significant loci to the strucuture
      if (!is.null(link2gene) & any(!is.null(threshold), is.null(threshold))) {
        link_association_structure <- transcript_association[transcript_association$Marker %in% 
                                                               link2gene$rs, ]
        if (dim(link_association_structure)[1] == 0) {
          print("no matched locis, will not draw linking line")
          link_asso_gene <- list(NULL)
          threshold_line <- list(NULL)
        }
        if (dim(link_association_structure)[1] > 0) {
          link_number <- dim(link_association_structure)[1]
          link_asso_gene <- list(geom_segment(data = link_association_structure, 
                                              aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, 
                                                                                 link_number), yend = -log10(p) * fold), linetype = "longdash", 
                                              colour = "red"))
        }
      }
      if (is.null(link2gene) & is.null(threshold)) {
        print("threshold acquired")
        link_asso_gene <- list(NULL)
      }
      if (is.null(link2gene) & !is.null(threshold)) {
        link_association_structure <- transcript_association[-log10(transcript_association$p) >= 
                                                               threshold, ]
        link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), 
        ]
        if (dim(link_association_structure)[1] == 0) {
          print("no -log10(p) pass the threshold, will not draw link")
          link_asso_gene <- list(NULL)
          threshold_line <- list(NULL)
        }
        if (dim(link_association_structure)[1] > 0) {
          link_association_structure <- transcript_association[-log10(transcript_association$p) >= 
                                                                 threshold, ]
          link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), 
          ]
          link_number <- dim(link_association_structure)[1]
          link_asso_gene <- list(geom_segment(data = link_association_structure, 
                                              aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, 
                                                                                 link_number), yend = -log10(p) * fold), linetype = "longdash", 
                                              colour = "red"))
        }
      }
      # add linking line to link gene structure and LD matrix
      if (isTRUE(triangleLD)) {
        if (is.null(link2gene) & is.null(link2LD)) {
          link_association_structure <- transcript_association[-log10(transcript_association$p) >= 
                                                                 threshold, ]
          link_association_structure <- link_association_structure[!duplicated(link_association_structure$p), 
          ]
          link_number <- dim(link_association_structure)[1]
          link_asso_gene <- list(geom_segment(data = link_association_structure, 
                                              aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, 
                                                                                 link_number), yend = -log10(p) * fold), linetype = "longdash", 
                                              colour = "red"))
          marker_axis_LD_x <- transcript_min + (seq(1:marker_number) - 
                                                  1) * 2 * distance
          marker_axis_genic_x <- hapmap_ld$pos
          marker_axis_LD_y <- rep(-5 * max(pvalue_range) * fold/30, marker_number)
          marker_axis_genic_y <- rep(-max(pvalue_range) * fold/30, marker_number)
          link_ld_data <- data.frame(x = marker_axis_LD_x, xend = marker_axis_genic_x, 
                                     y = marker_axis_LD_y, yend = marker_axis_genic_y)
          link_ld_data <- link_ld_data[link_ld_data$xend %in% link_association_structure$Site, 
          ]
          link_LD_genic_structure <- geom_segment(data = link_ld_data, 
                                                  aes(x = x, xend = xend, y = y, yend = yend), colour = "red", 
                                                  linetype = "longdash")
        }
        if (!is.null(link2gene) & !is.null(link2LD)) {
          link_association_structure <- transcript_association[transcript_association$Marker %in% 
                                                                 link2LD$rs, ]
          link_number <- dim(link_association_structure)[1]
          link_asso_gene <- list(geom_segment(data = link_association_structure, 
                                              aes(x = Site, xend = Site, y = rep(-max(pvalue_range) * fold/30, 
                                                                                 link_number), yend = -log10(p) * fold), linetype = "longdash", 
                                              colour = "red"))
          marker_axis_LD_x <- transcript_min + (seq(1:marker_number) - 
                                                  1) * 2 * distance
          marker_axis_genic_x <- hapmap_ld$pos
          marker_axis_LD_y <- rep(-4.5 * max(pvalue_range) * fold/30, marker_number)
          marker_axis_genic_y <- rep(-max(pvalue_range) * fold/30, marker_number)
          link_ld_data <- data.frame(x = marker_axis_LD_x, xend = marker_axis_genic_x, 
                                     y = marker_axis_LD_y, yend = marker_axis_genic_y)
          link_ld_data <- link_ld_data[link_ld_data$xend %in% link_association_structure$Site, 
          ]
          link_LD_genic_structure <- geom_segment(data = link_ld_data, 
                                                  aes(x = x, xend = xend, y = y, yend = yend), colour = "red", 
                                                  linetype = "longdash")
        }
      }
      if (!is.null(link2gene) & is.null(link2LD)) {
        link_LD_genic_structure <- list(NULL)
      }
      if (is.null(link2gene) & !is.null(link2LD)) {
        link_LD_genic_structure <- list(NULL)
      }
      if (!isTRUE(triangleLD)) {
        link_LD_genic_structure <- list(NULL)
      }
      y_axis_text <- ifelse(unique(transcript_corrdination$V7) == "-", list(geom_text(aes(x = transcript_max + 
                                                                                            (transcript_max - transcript_min)/6.5, y = mean(pvalue_range) * 
                                                                                            fold), label = "atop(-log[10]*italic(P)[observed])", parse = TRUE, 
                                                                                      angle = 90)), list(geom_text(aes(x = transcript_min - (transcript_max - 
                                                                                                                                               transcript_min)/6.5, y = mean(pvalue_range) * fold), label = "atop(-log[10]*italic(P)[observed])", 
                                                                                                                   parse = TRUE, angle = 90)))
      if (isTRUE(triangleLD)) {
        xtext <- list(geom_text(aes(x = (transcript_max + transcript_min)/2, 
                                    y = min(poly_data$y) - 10 * distance, label = transcript)))
      } else {
        xtext <- list(geom_text(aes(x = (transcript_max + transcript_min)/2, 
                                    y = -(transcript_max - transcript_min)/10, label = transcript)))
      }
      # add shape for the points of leadsnp
      leadsnp2highlight <- transcript_association[transcript_association$Marker == 
                                                    leadsnp, ]
      leadsnp2highlight_list <- list(geom_point(data = leadsnp2highlight, 
                                                aes(x = Site, y = -log10(p) * fold), shape = leadsnp_shape, colour = leadsnp_colour, 
                                                fill = leadsnp_fill, size = leadsnp_size))
      # change the shape,size,colour, and fill for highlighted marker
      if (is.null(marker2highlight)) {
        marker2highlight_list = list(NULL)
      } else {
        marker2highlight <- merge(marker2highlight, transcript_association, 
                                  by.x = "rs", by.y = "Marker")
        # marker2highlight_list = list(geom_point(data=marker2highlight, aes(Site,
        # -log10(p) * fold, shape=factor(shape), colour=factor(colour),
        # fill=factor(fill), size=factor(size))))
        
        marker2highlight_list = list(annotate("point", x = marker2highlight$Site, 
                                              y = -log10(marker2highlight$p) * fold, shape = marker2highlight$shape, 
                                              colour = marker2highlight$colour, size = marker2highlight$size, 
                                              fill = marker2highlight$fill))
      }
      if (!is.null(marker2label)) {
        marker2label <- merge(marker2label, transcript_association, by.x = "rs", 
                              by.y = "Marker")
        # marker2label_list <-
        # list(annotate('text',x=marker2label$Site,y=-log10(marker2label$p) *
        # fold,label=marker2label$rs,angle=marker2label_angle))
        marker2label_list <- list(geom_text_repel(aes(x = marker2label$Site, 
                                                      y = -log10(marker2label$p) * fold, label = marker2label$rs), 
                                                  angle = marker2label_angle, size = marker2label_size))
      } else {
        marker2label_list <- list(NULL)
      }
      # plot the reduced point if highlighted marker exited
      if (is.null(marker2highlight)) {
        transcript_association = transcript_association
      } else {
        transcript_association = transcript_association[!(transcript_association$Marker %in% 
                                                            marker2highlight$rs), ]
      }
      plot <- ggplot() + geom_point(data = transcript_association, aes(Site, 
                                                                       -log10(p) * fold), colour = "black",size=upperpointsize) + ld_leadsnp_colour + transcript_intron_structure + 
        transcript_structure_exon_list + transcript_structure_utr_list + 
        transcript_structure_cds_list + link_asso_gene + link_LD_genic_structure + 
        scale_x + scale_y_line + scale_y_ticks + scale_y_text + threshold_line + 
        bottom_trianglLD + y_axis_text + xtext + leadsnp2highlight_list + 
        marker2label_list + marker2highlight_list + theme_bw() + theme(legend.key = element_rect(colour = "black"), 
                                                                       axis.title.x = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), 
                                                                       panel.grid = element_blank(), axis.text = element_blank(), axis.title.y = element_blank(), 
                                                                       text = element_text(size = 15, face = "bold"))
      return(plot)
    }
  }
}
