## date: 2021-01-12
## author: alice yue
## input: 2D matrices (cell x marker) + their clr (cell x cell population true/false matrix)
## ouput: gigaSOM clustering

# i start my docker from the designated folder already, cd there if you don't
# sudo docker run -it --rm -v "$PWD":/mnt/FCS_local3/backup/'Brinkman group'/current/Alice/flowMagic_data -w /mnt/FCS_local3/backup/'Brinkman group'/current/Alice/flowMagic_data julia

# install and load package
import Pkg
Pkg.add("GigaSOM")
Pkg.add("CSVFiles")
Pkg.add("Glob")
Pkg.add("DataFrames")
# Pkg.test("GigaSOM"); # test installation
using GigaSOM
using CSVFiles, DataFrames
using Glob # load package

# enable @debug messages
using Logging
global_logger(ConsoleLogger(stderr, Logging.Debug))

# get files
cd("/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data")

datasets = readdir("data/2D/x")
for dataset in datasets
  scats = readdir("data/2D/x/$dataset")
  for scat in scats
    csvfiles = Glob.glob("*.csv.gz","data/2D/x/$dataset/$scat")
    for csvfile in csvfiles
      exprs = DataFrame(load(File(format"CSV", csvfile)))
      # clr = DataFrame(load(File(format"CSV", replace(csvfile,"/x/"=>"/y/"))))
      # ncpop = ncol(clr)
      # ncpop = names(clr)[end]=="other" ? ncpop-1 : ncpop

      som = GigaSOM.initGigaSOM(exprs, 2, 2)    # 4 clusters! at most 4 clusters in a scatterplot
      som = GigaSOM.trainGigaSOM(som, exprs)      # SOM training
      clusters = GigaSOM.mapToGigaSOM(som, exprs) # extraction of per-cell cluster IDs
      # e = GigaSOM.embedGigaSOM(som, exprs)        # EmbedSOM projection to 2D

      save(File(format"CSV", replace(csvfile,"data/2D/x"=>"results/2D/GigaSOM_clusters")), clusters)
    end
  end
end

datasets = readdir("data/nD/x")
for dataset in datasets
  csvfiles = Glob.glob("*.csv.gz","data/nD/x/$dataset")
  for csvfile in csvfiles
    exprs = DataFrame(load(File(format"CSV", csvfile)))
    clr = DataFrame(load(File(format"CSV", replace(csvfile,"/x/"=>"/y/"))))
    ncpop = ncol(clr)
    ncpop = names(clr)[end]=="other" ? ncpop-1 : ncpop

    ncpop1 = Int(ceil(sqrt(ncpop)))

    som = GigaSOM.initGigaSOM(exprs, ncpop1, ncpop1)
    som = GigaSOM.trainGigaSOM(som, exprs)      # SOM training
    clusters = GigaSOM.mapToGigaSOM(som, exprs) # extraction of per-cell cluster IDs
    # e = GigaSOM.embedGigaSOM(som, exprs)        # EmbedSOM projection to 2D

    save(File(format"CSV", replace(csvfile,"data/nD/x"=>"results/nD/GigaSOM_clusters")), clusters)
  end
end

# going to use a matrix, not fcs
# params, fcsmatrix = loadFCS("Levine_13dim.fcs")  # load the FCS file
# exprs = fcsmatrix[:,1:13]  # extract only the data columns with expression values
