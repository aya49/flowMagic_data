# flowMagic_data

## a few notes on the scripts
- 00_generate_data_<data set>:
  - <hipc> data set has 2 panels where until "Live cells" the gating is the same.
  - straight thresholds for <hipc> are usually min (so - is <, + is >=); if not, cell populations are already adjusted accordingly

## deepCyTOF

in order to get the original code to work with minimal modifications, we used:
- scikit-learn 0.20.0
- keras 1.2.2
- TensorFlow 1.0.0
- numpy 1.19

## folder structure

```{bash}
├─── src (source code)
├──+ data (raw data)
|  ├─+ 2D (2D scatterplots)
|  | ├── x/<data_set>/<scatterplot>/<FCS name>.csv.gz (cell x 2 markers) cleaned transformed FI value matrices
|  | ├── y/<data_set>/<scatterplot>/<FCS name>.csv.gz (cell x cell population) 0/1 label matrices
|  | └── thresholds/<data_set>/<FCS name>.Rdata (list of scatterplots > markers > threshold value)
|  └─+ nD (full nD FCS file)
|    ├── x/<data_set>/<FCS name>.csv.gz (cell x markers) cleaned transformed FI value matrices
|    ├── y/<data_set>/<FCS name>.csv.gz (cell x leaf cell population) 0/1 label matrices
|    └── thresholds/<data_set>/<FCS name>.Rdata (list of scatterplots > list of markers for each scatterplot > threshold value)
├──+ results
|  ├─+ 2D (*input data for flowMagic: 400 x 400)
|  | ├── deepCyTOF_labels/<data_set>/<scatterplot>/<FCS name>.csv (cell vector) int cell population labels corresponding to y column (0 means none)
|  | ├── deepCyTOF_models/<data_set>/<scatterplot>/cellClassifier.h5 (model files)
|  | ├── flowLearn_thresholds/<data_set>/<scatterplot>/<k>.Rdata (list of markers for each scatterplot > threshold value)
|  | ├── flowLearn_plots/<data_set>/<scatterplot>/<k>/<FCS name>(_train).png ("_train" if the sample is a reference sample for at least 1 threshold)
|  | ├── GigaSOM_clusters/<data_set>/<scatterplot>/<FCS name>.csv.gz (cell x 1) int cluster labels (2x2=4 clusters)
|  | ├── GigaSOM_plots/<data_set>/<scatterplot>/<FCS name>.png (2 plots, actual convex hull vs clustered)
|  | ├─* x_2Ddensity/<data_set>/<scatterplot>/<FCS name>.csv.gz (smoothed 2D Gaussian kernel density)
|  | ├─* x_2Dscatter/<data_set>/<scatterplot>/<FCS name>.csv.gz (1/0 whether or not there is a cell on pixel)
|  | ├─* x_2Dncells/<data_set>/<scatterplot>/<FCS name>.csv.gz (number of cells in each pixel)
|  | ├─* y_2D/<data_set>/<scatterplot>/<FCS name>.csv.gz (string cell population label of each pixel)
|  | └── y_vector/<data_set>/<scatterplot>/<FCS name>.csv.gz (cell vector) string cell population label of each cell, if needed
|  └─+ nD
|    ├── deepCyTOF_labels/<data_set>/<FCS name>.csv (cell x markers) (cell vector) int cell population labels corresponding to y column (0 means none)
|    ├── deepCyTOF_models/<data_set>/cellClassifier.h5 (model files)
|    └── GigaSOM_clusters/<data_set>/<FCS name>.csv.gz (cell x 1) int cluster labels (hxh=g clusters >= actual number of cell populations)
└──+ scores
   ├─+ 2D
   | ├── deepCyTOF/<data_set>/<scatterplot>/<FCS name>.csv.gz (F1 score for all cpops)
   | ├── flowLearn/<data_set>/<scatterplot>/<k>.csv.gz (F1 scores for all FCS files > all cpops)
   | ├── GigaSOM/<data_set>/<scatterplot>/<FCS name>.csv.gz (F1 score for all cpops)
   | ├── 
   | ├── 
   | ├── 
   | ├── 
   | └── 
   └─+ nD
     ├── deepCyTOF/<data_set>/<FCS name>.csv.gz (F1 score for all leaf cpops)
     ├── 
     └── 
```
