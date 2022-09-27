library(plyr)

#compute into groups by promoter methylation (By eye according to WGBS heatmap)
Stem_Cells <- WGBS_data[,1:10]
Stem_Cells <- rowMeans(Stem_Cells, na.rm = TRUE)

Organs <- WGBS_data[,c('Gastric',
                        'Left Ventricle',
                        'Psoas Muscle',
                        'Right Atrium',
                        'Right Ventricle',
                        'Sigmoid Colon',
                        'Small Intestine',
                        'Esophagus')]
Organs <- rowMeans(Organs, na.rm = TRUE)
Neuro_cells <- WGBS_data[,c('Neurosphere Cultured Cells Cortex Derived'
                            ,'Neurosphere Cultured Cells Ganglionic Eminence Derived'
                            ,'Brain Germinal Matrix'
                            ,'Brain Hippocampus Middle')]
Neuro_cells <- rowMeans(Neuro_cells, na.rm = TRUE)

Ovary_Aorta_Psoas <- WGBS_data[,c('Ovary','Aorta','Psoas Muscle')]
Ovary_Aorta_Psoas <- rowMeans(Ovary_Aorta_Psoas, na.rm = TRUE)

Fetal_Intestine_Large <- WGBS_data[,'Fetal Intestine Large']

Fetal_Intestine_Small <- WGBS_data[,'Fetal Intestine Small']

Liver_Pancreas <- WGBS_data[,c('Adult Liver', 'Pancreas')]
Liver_Pancreas <- rowMeans(Liver_Pancreas, na.rm = TRUE)

Foreskin_Keratinocyte <- WGBS_data[,'Penis Foreskin Keratinocyte Primary Cells skin03']

Immune_system <- WGBS_data[,c('Thymus','Spleen','Mobilized CD34 Primary Cells Female')]
Immune_system <- rowMeans(Immune_system, na.rm = TRUE)

Grouped_WGBS <- data.frame(Stem_Cells, Organs, Neuro_cells, Ovary_Aorta_Psoas, Fetal_Intestine_Large, Fetal_Intestine_Small, Liver_Pancreas, Foreskin_Keratinocyte, Immune_system)
write.csv(Grouped_WGBS, 'Grouped_WGBS.csv')

#t-test

combos <- combn(ncol(Grouped_WGBS),2)
adply(combos, 2, function(x) {
  test <- t.test(Grouped_WGBS[, x[1]], Grouped_WGBS[, x[2]])
  
  out <- data.frame("var1" = colnames(Grouped_WGBS)[x[1]]
                    , "var2" = colnames(Grouped_WGBS[x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
})
