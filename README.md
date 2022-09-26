# flowMagic_data

document: https://www.overleaf.com/4958349237rghgfgvdmrxq

## a few notes on the scripts
- 00_generate_data_<data set>:
  - <hipc> data set has 2 panels where until "Live cells" the gating is the same.
  - straight thresholds for <hipc> are usually min (so - is <, + is >=); if not, cell populations are already adjusted accordingly

## deepCyTOF

we used the original package versions used in the deepCyTOF paper to run deepCyTOF: 
- scikit-learn 0.20.0
- keras 1.2.2
- TensorFlow 1.0.0
- numpy 1.19

scripts to run deepCyTOF on current versions of packages are available in the [Util](Util) folder; these scripts end in `_.py`.

## folder structure

2D data sets: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7108752.svg)](https://doi.org/10.5281/zenodo.7108752)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7108766.svg)](https://doi.org/10.5281/zenodo.7108766)

```{bash}
<dsfkm> = <data_set>/<scatterplot>/<FCS name>/<number of training samples used>/<marker>
* [400 x 400] in/output data for flowMagic

├─── src   (source code)
├──+ raw   (raw data)
|  ├─+ 2D   (2D scatterplots)
|  | ├── x/<dsf>.csv.gz   [cell x 2 markers] (cleaned transformed FI value matrices)
|  | ├── y/<dsf>.csv.gz   [cell x cell population] (0/1 label matrices)
|  | ├── filters/<dsf>.Rdata   [list of cell populations > convex hull: points x 2 markers matrix (first == last point)] (used for labelling cell populations on plots)
|  | └── thresholds/<dsf>.Rdata   [list of scatterplots > markers > threshold numeric value]
|  ├─+ nD   (full nD FCS file)
|  | ├── x/<df>.csv.gz   [cell x markers] (cleaned transformed FI value matrices)
|  | └── y/<df>.csv.gz   [cell x leaf cell population] (0/1 label matrices)
|  ├─+ scatterplots/<dsf>.png   (flowDensity gating scatterplots)
|  ├─+ HIPCbcell.pdf   (manual gating hierarchy tree for data set HIPCbcell)
|  ├─+ HIPCmyeloid.pdf   (manual gating hierarchy tree for data set HIPCmyeloid)
|  └─+ pregnancy_dates.png   (density kernel estimate for all threshold gates across all samples in the pregnancy data set)
├──+ data    (preprocessed data)
|  ├─+ 2D
|  | ├─* x_2Ddensity/<dsf>.csv.gz   (smoothed 2D Gaussian kernel density)
|  | ├─* x_2Dscatter/<dsf>.csv.gz   (1/0 whether or not there is a cell on pixel)
|  | ├─* x_2Dncells/<dsf>.csv.gz   (number of cells in each pixel)
|  | ├── x_2Ddensity_euclidean/<ds>.Rdata   (pair-wise sample distance object)
|  | ├── x_2Ddensity_euclidean_rankkmed/<ds>/<k>/<reference FCS name>.csv.gz   (vector of FCS sample names in the said cluster)
|  | ├─* y_2D/<dsf>.csv.gz   (numeric label matrix; 0 means other)
|  | ├─* y_2Dncells/<dsf>.csv.gz   (numeric count matrix indicating number of cells on each pixel)
|  | └── y_vector/<data_set>/<scatterplot>/<FCS name>.csv.gz   (cell vector) string cell population label of each cell, if needed
|  └─+ nD
|    └── Rtsne/<df>.csv.gz   [cell x 2] (2D Rtsne dimensionality reduced nD FCS files for plotting)
├──+ results
|  ├─+ 2D
|  | ├─+ deepCyTOF_F1/<kds>
|  | | ├── F1.csv   (vector of cell population count weighted F1 score (calculated by deepCyTOF) for each FCS sample)
|  | | └── fnames.csv   (FCS names for F1.csv^)
|  | ├── deepCyTOF_labels/<kdsf>.csv   (cell vector of int cell population labels corresponding to y column; 0 means other/none)
|  | ├── deepCyTOF_models/<kds>/cellClassifier.h5   (model files)
|  | ├── flowLearn_labels/<dskf>.csv.gz   [cell x cell population] (TRUE/FALSE cell population label matrices)
|  | ├── flowLearn_thresholds/<dsmk>.Rdata   (list of markers for each scatterplot > threshold value)
|  | ├── flowLearn_plots/<dskf>.png   (these plots are compiled in plots/2D/all/flowLearn)
|  | ├── gigaSOM_clusters/<dsf>.csv.gz   (cell vector of int cluster labels (2x2=4 clusters))
|  | └── gigaSOM_lbls/<dsf>.csv.gz   [cell x cell population] (TRUE/FALSE cell population label matrices) !!!will change this to labels (not lbls)
|  └─+ nD
|    ├─+ deepCyTOF_F1/<d>
|    | ├── F1.csv   (vector of cell population count weighted F1 score (calculated by deepCyTOF) for each FCS sample)
|    | └── fnames.csv   (FCS names for F1.csv^)
|    ├── deepCyTOF_labels/<df>.csv   (cell vector of int cell population labels corresponding to y column (0 means other/none)
|    ├── deepCyTOF_models/<d>/cellClassifier.h5   (model files)
|    ├── gigaSOM_clusters/<df>.csv.gz   (cell vector of int cluster labels (hxh=g clusters >= actual number of cell populations))
|    └── gigaSOM_labels/<df>.csv.gz   [cell x cell population] (TRUE/FALSE cell population label matrices)
├──+ scores   (F1 score for each cpop and FCS sample)
|  ├─+ 2D   [. x method, dataset, scatterplot, cpop, train_no, fcs, train, precision, recall, f1, true_proportion, predicted_proportion, true_size, predicted_size]
|  | ├── deepCyTOF/<kds>.csv.gz
|  | ├── flowLearn/<ds>.csv.gz
|  | └── gigaSOM/<ds>.csv.gz
|  ├─+ nD
|  | ├── deepCyTOF/<kd>.csv.gz
|  | ├── flowLearn/<d>.csv.gz
|  | └── gigaSOM/<d>.csv.gz
|  ├── 2D.csv.gz   (compiled matrix from 2D folder)
|  └── nD.csv.gz   (compiled matrix from nD folder)
└──+ plots
   ├─+ 2D 
   | ├─+ all   (compiled scatterplots; coloured points are results, convex hulls are original density gates)
   | | ├── original/<ds>/<1,2,...>.png
   | | ├── deepCyTOF/<ds>/<1,2,...>.png
   | | ├── gigaSOM/<ds>/<1,2,...>.png
   | | └── flowLearn/<ds>/<1,2,...>.png
   | └─+ scores    (plots relating to F1 scores etc.)
   └─+ nD 
```
