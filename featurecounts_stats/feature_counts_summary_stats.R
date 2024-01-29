dd <- read.csv('counts.txt.summary', sep = '\t')


#get into R format
dd_transpose <- t(dd)
colnames(dd_transpose) <- dd_transpose[1, ]
dd_transpose <- dd_transpose[-1,]
dd <- dd_transpose
dd <- as.data.frame(dd)

#convert character to numeric type
dd <- as.data.frame(lapply(dd, as.numeric))


#get column for total counts 
dd$Total <- dd$Assigned + dd$Unassigned_Unmapped + dd$Unassigned_MultiMapping + dd$Unassigned_NoFeatures + dd$Unassigned_Ambiguity

#how many are assigned?
dd$prop.assigned <- (dd$Assigned / dd$Total)*100

print(mean(dd$prop.assigned))
sd(dd$prop.assigned)

#how many are ambiguous?
dd$prop.ambiguous <- (dd$Unassigned_Ambiguity /dd$Total)*100
mean(dd$prop.ambiguous)
sd(dd$prop.ambiguous)
