
library(GEOquery)

GSE38959 <- getGEO('GSE38959', destdir=".")

info = GSE38959[["GSE38959_series_matrix.txt.gz"]]@featureData@data

x = data.frame(GSE38959[["GSE38959_series_matrix.txt.gz"]]@assayData[["exprs"]])

# eliminar últimas 17 colunas que não são de TNBC
x = x[,-31:-47]

# nome dos genes
gene = info$ENSEMBL_ID #matlab

gene = info$GENE_SYMBOL#python

x = cbind(gene, x)

# apagar as linhas que não tenham o nome alternativo do gene
x = x[x$gene != "",]

# apaga genes duplicados
x = x[!duplicated(x[,c('gene')]),]



# usar apenas n linhas:
n = dim(x)[1]
x = x[1:n,]


row.names(x)=x$gene

x[c('gene')] <- list(NULL)

#Normalização
mean=rowMeans(x)
expr=5*log(1+(x/mean))

express_T=t(expr) #transposto do dataframe

write.csv(express_T, 'D:\\universidade\\projeto\\r\\expression.csv')
write.csv(expr, 'D:\\universidade\\projeto\\r\\expression2.csv')