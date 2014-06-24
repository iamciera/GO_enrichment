library(goseq)
library(GO.db)

pat <- read.table("tf2cmbr_wtcmbr_DE1_sigonly_GOmatrix.txt",header=TRUE)
cate <- read.table("melted.GOTable.txt",header=TRUE)


genes = as.integer(pat$all)
names(genes) = pat$itag
table(genes)
length(genes)

pwf = nullp(genes,bias.data=pat$length)


GO.wall = goseq(pwf,gene2cat = cate)
GO.wall

enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05]

enriched.GO

my.GO <- as.character(enriched.GO)
my.GO.table <- Term(my.GO)
my.GO.table
t <- as.matrix(my.GO.table)

write.table(t, file="tf2cmbr_wtcmbr_DE1_sigonly_ALL_GOmerge.txt")


