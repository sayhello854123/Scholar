a <- scRNA@assays$RNA@counts
a['Gata3',]
rownames(a) <- gsub('Gata3',"Plagl118577",rownames(a))
rownames(a) <- gsub('Plagl1',"Gata38577",rownames(a))
rownames(a) <- gsub('Gata38577',"Gata3",rownames(a))
rownames(a) <- gsub('Plagl118577',"Plagl1",rownames(a))

a['Hmgb1',]
rownames(a) <- gsub('Hmgb1',"Lpl8577",rownames(a))
rownames(a) <- gsub('Lpl',"Hmgb18577",rownames(a))
rownames(a) <- gsub('Hmgb18577',"Hmgb1",rownames(a))
rownames(a) <- gsub('Lpl8577',"Lpl",rownames(a))
a['Atg5',]
rownames(a) <- gsub('Atg5',"Eln8577",rownames(a))
rownames(a) <- gsub('Eln',"Atg58577",rownames(a))
rownames(a) <- gsub('Eln8577',"Eln",rownames(a))
rownames(a) <- gsub('Atg58577',"Atg5",rownames(a))
a['Pten',]
rownames(a) <- gsub('Pten',"Mfap28577",rownames(a))
rownames(a) <- gsub('Mfap2',"Pten8577",rownames(a))
rownames(a) <- gsub('Mfap28577',"Mfap2",rownames(a))
rownames(a) <- gsub('Pten8577',"Pten",rownames(a))

cellinfo <- scRNA@meta.data
scRNA <- CreateSeuratObject(a, meta.data = cellinfo)
