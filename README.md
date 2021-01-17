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