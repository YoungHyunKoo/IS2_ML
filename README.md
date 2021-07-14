# IS2_ML
ICESat-2 sea ice processing for surface classification and freeboard calculation

- 01_find_overlapped_data.ipynb: Find overlapped ICESat-2 and Sentinel-2 data, including download ICESat-2 data.
- 02_ATL07_processing.ipynb: Process downloaded ATL07 files and save them as shapefiles and csvfiles.
- visualize_images.py: Display ATL07 shapefiles and coincident Sentinel-2 images via QGIS software. Label ATL07 data using ground-truth Sentinel-2 images
- 03_ATL07_train_ANN.ipynb: Train and test Machine Learning models (e.g. DNN and LSTM) based on the labelled train data. 
