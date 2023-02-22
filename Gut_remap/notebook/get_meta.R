library(Seurat)
path = '/domino/edv/id-td-virology/Public_dataset/2022_Gut/aggr/analysis/20220817_subset_full_data.B_seurat_object_afterresolution_optimization.rds'
seu <- readRDS(file = path)
View(seu@meta.data)

savepath = '/domino/edv/id-td-virology/Zhiyuan/public/Gut_remap/processed_data/B_meta.tsv'
write.table(seu@meta.data, file = savepath, quote = FALSE, sep = '\t')





