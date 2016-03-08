
function is_cosmicgene(val)
    length(val) == 0 ? false : true
end

function is_highfreqgene(gene)
    in(gene, Set{ASCIIString}(["TP53", "EGFR", "KRAS","CDKN2A","KMT2C","STK11","NF1","KMT2D","GRIN2A","SMARCA4",#]))
                               "ZNF521","RB1","SETBP1","KDR","ATM","NFE2L2","CREBBP","ALK","PIK3CA","TRRAP"]))
end

function readlung()
    data = readcsv("/haplox/ref/guolin_probe.csv",ASCIIString)
    lung_genes = filter(x->length(x)==0?false:true, map(strip,data[2:end,1]))
    tuhao_genes = filter(x->length(x)==0?false:true, map(strip, data[2:end,end]))
    @show length(lung_genes)
    @show length(tuhao_genes)
    @show tuhao_genes
    @show "TP53" in tuhao_genes
    lung = filter(is_highfreqgene,lung_genes)
    tuhao = filter(is_highfreqgene, tuhao_genes)
    @show tuhao
    @show length(lung)
    @show length(tuhao)
end

function in_cosmic(data,csv)
    raw_genes = Set{ASCIIString}(data[2:end,6])
#    info("There are $(length(raw_genes)) raw genes in $csv")
    cosmic_id = data[:, end]
    
    inds_cosmic = map(is_cosmicgene, cosmic_id)
    gene = data[:,6]
    
    gene_cosmic = gene[inds_cosmic,:]
    cosmic_genes = Set{ASCIIString}(gene_cosmic[2:end])
#    info("There are $(length(cosmic_genes)) cosmic genes in $csv")
    #@show gene_cosmic[1:5]
    #break
#    info(repeat("=",70))
    
    cosmic_genes
end

function get_goodfeature(data_dir)
    csvs = map(x->joinpath(data_dir,x), readdir(data_dir))
    csvs = filter(x->endswith(x,".csv"), csvs)
    cosmic_genes = Set{ASCIIString}()
 #  features = Array{Int64,2}(0,)
    i =  1
    for csv in csvs
        data = readcsv(csv, ASCIIString)
        cosmic_genes = in_cosmic(data,csv)
        genes = filter(is_highfreqgene, cosmic_genes)
        @show genes
        info("There are $(length(genes)) hight frequent cosmic genes in $csv")
        i % 2 == 0 && info(repeat("=",60))
        i+=1
    end
end

using TSne

function cluster_tsne(data,label)
    Xcenter = data - mean(data)
    Xstd = std(data)
    X = Xcenter / Xstd
    Y = tsne(X, 2, 50, 1000, 20.0)
    
    labels = ASCIIString[string(l) for l in label]
    theplot = plot(x=Y[:,1], y= Y[:,2], color=labels)

    draw(PDF("tnse_plot.pdf", 4inch, 3inch), theplot)
end

using MultivariateStats

function cluster_pca(data)
    M = fit(PCA,data, maxoutdim=2)
    transform(M, data)
end


