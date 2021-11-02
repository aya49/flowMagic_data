cd /home/aya43/flowMagic_data/data

wget http://vault.sfu.ca/index.php/s/1zDFwAbR39MFAP5/download -O y_vector_.zip
unzip y_vector_.zip
rm y_vector_.zip

mkdir x_2Ddiscrete
cd x_2Ddiscrete

wget http://vault.sfu.ca/index.php/s/nK9p7f3YGAS8QAh/download -O sangerP2.zip
unzip sangerP2.zip
rm sangerP2.zip

wget http://vault.sfu.ca/index.php/s/32V5vXMpsMDWJjX/download -O HIPCmyeloid_pregnancy.zip
unzip HIPCmyeloid_pregnancy.zip
rm HIPCmyeloid_pregnancy.zip

wget http://vault.sfu.ca/index.php/s/25ZNncQWeLYagql/download -O HIPCbcell.zip
unzip HIPCbcell.zip
rm HIPCbcell.zip

cd ..

wget http://vault.sfu.ca/index.php/s/4vcFP82sIR6S53H/download -O y_2D.zip
unzip y_2D.zip
rm y_2D.zip

wget http://vault.sfu.ca/index.php/s/ursZXAFtxEG9kVI/download -O x_2Ddenscat.zip
unzip x_2Ddenscat.zip
rm x_2Ddenscat.zip

wget http://vault.sfu.ca/index.php/s/MzA6acXFQYNCO91/download -O x_2Dcontour.zip
unzip x_2Dcontour.zip
rm x_2Dcontour.zip

wget http://vault.sfu.ca/index.php/s/09lBvEI1w1ao5KD/download -O x_2Ddenscat_euclidean_rankkmed.zip
unzip x_2Ddenscat_euclidean_rankkmed.zip
rm x_2Ddenscat_euclidean_rankkmed.zip

cd /home/aya43/flowMagic_data

mkdir raw
cd raw
mkdir 2D
cd 2D
wget http://vault.sfu.ca/index.php/s/dVnVhG1nxg2x2bN/download -O y.zip
unzip y.zip
rm y.zip




# miniconda
conda install -c anaconda numpy
conda install -c conda-forge pandas
pip install mmcv-full==1.3.14 -f https://download.openmmlab.com/mmcv/dist/cu101/torch1.7.0/index.html
pip install mmsegmentation

pip install argparse
pip install compress_pickle
pip install https://github.com/ufoym/imbalanced-dataset-sampler/archive/master.zip
pip install tensorboard_logger
pip install GPUtil

