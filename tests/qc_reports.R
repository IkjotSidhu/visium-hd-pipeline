# ---- QC reporting: weighted means, independent groups, empty entities ----
source('R/utils.R')
.libPaths(c(default_config$r_library,.libPaths()))
suppressPackageStartupMessages(library(tidyverse))
cfg <- default_config; cfg$bin_size <- '8'; cfg$res_label <- '008um'
cfg$qc$min_counts <- 50; cfg$qc$min_features <- 5; cfg$qc$max_mt_pct <- 50
cfg$qc$report$sample_col <- 'sample'; cfg$qc$report$group_cols <- c('condition','donor')
sheet <- tibble(replicate_id=c('R1','R2','R3','R4'), sample=c('S1','S1','S2','S2'),
 library_id=c('L1','L1','L2','L2'), donor=c('D1','D1','D2','D2'),condition=c('A','A','B','B'))
data <- sheet[c(1,1,2,3,3),] |> mutate(UMI=c(20,60,100,80,1000),genes=c(10,30,50,100,1),
 percent.mt=c(5,10,20,80,5),pass=c(FALSE,TRUE,TRUE,FALSE,FALSE))
reports <- build_qc_reports(data,sheet,cfg,'UMI','genes')
a <- reports$group$metrics |> filter(group_variable=='condition',group=='A',stage=='before_qc')
stopifnot(a$mean_umi == 60, a$n_units == 3) # Not the unweighted replicate mean 70.
after <- reports$sample$summary |> filter(sample_id=='S1')
stopifnot(after$mean_umi == 80,after$mean_genes == 40,after$total_umi == 160,
          after$n_before == 3,after$n_after == 2,abs(after$pct_kept-200/3)<1e-9)
empty <- reports$sample$summary |> filter(sample_id=='S2')
stopifnot(empty$n_after==0,empty$total_umi==0,is.na(empty$mean_umi),empty$pct_kept==0)
absent <- reports$replicate$summary |> filter(replicate_id=='R4')
stopifnot(absent$n_before==0,absent$n_after==0,is.na(absent$pct_kept))
b <- reports$group$retention |> filter(group_variable=='condition',group=='B')
stopifnot(b$fail_mt==1,b$fail_genes==1,b$fail_any==2)
# Overlapping flags count one removed unit only.
data$genes[4] <- 1
reports2 <- build_qc_reports(data,sheet,cfg,'UMI','genes')
b <- reports2$group$retention |> filter(group_variable=='condition',group=='B')
stopifnot(b$fail_mt+b$fail_genes==3,b$fail_any==2)
# Missing group metadata remains a visible group, not dropped observations.
sheet$condition[1] <- NA_character_; data$condition[1:2] <- NA_character_
missing <- build_qc_reports(data,sheet,cfg,'UMI','genes')
stopifnot(any(is.na(missing$group$summary$group)))
root <- tempfile('vhd-qc-'); dir.create(root)
cfg$step_dir <- root; cfg$project_name <- '<QC test>'; cfg$test_fraction <- 0.1
reports$test_fraction <- 0.1; dir.create(file.path(root,'tables'))
write_qc_reports(reports,cfg)
html <- paste(readLines(file.path(root,'qc_report.html')),collapse='\n')
stopifnot(all(c('qc_by_group.csv','qc_by_sample.csv','qc_by_replicate.csv') %in% list.files(file.path(root,'tables'))),
          grepl('Mean UMI/unit',html,fixed=TRUE),grepl('Test run: 10%',html,fixed=TRUE),
          grepl('&lt;QC test&gt;',html,fixed=TRUE),!grepl('<QC test>',html,fixed=TRUE))
unlink(root,recursive=TRUE)
cat('QC reporting checks passed\n')
