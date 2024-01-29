dd <- read.csv('mapping_data_rna_seq.csv', sep = ';')


give_sd_mean <- function(group, align) {
  mean_group <- mean(group[[align]])
  sd_group <- sd(group[[align]])
  
  print(paste0('Align: ', align, ', Mean: ', mean_group, ', Standard Deviation: ', sd_group))
}


# Define the four groups
blood_control <- dd[dd$group == 'Blood_Control', ]
blood_case <- dd[dd$group == 'Blood_Case', ]
lung_control <- dd[dd$group == 'Lung_Control', ]
lung_case <- dd[dd$group == 'Lung_Case', ]

alignments <- list('align_0', 'align_1', 'align_multi', 'overall')

for (element in alignments) {
  groups <- list(blood_control, blood_case, lung_control, lung_case)
  for (group in groups) {
    give_sd_mean(group, element)
  }
}


#what does it look like for all th samples combined'
mean(dd$align_0)
sd(dd$align_0)

mean(dd$align_1)
sd(dd$align_1)

mean(dd$align_multi)
sd(dd$align_multi)

mean(dd$overall)
sd(dd$overall)


#are there any outliers?
boxplot(dd$align_0)
boxplot(dd$align_1)
boxplot(dd$align_multi)
boxplot(dd$overall)






