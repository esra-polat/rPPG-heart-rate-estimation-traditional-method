> This project is the first phase of the rPPG Based Heart Rate Estimation Using Deep Learning, which is our graduation project. You can go to the full version of the project from this link.

#### Developers
* Esra POLAT - https://github.com/esra-polat
* Nur Deniz ÇAYLI - https://github.com/nurdenizcayli
* Minel SAYGISEVER - https://github.com/minelsaygisever

------------

Heart rate values were measured with UBFC dataset by adding face detection and skin segmentation to the 
- Chrominance-based Method (CHROM)
- Green - Vercruysse Method (GREEN)
- Independent Component Analysis Method (ICA)
- Plane-Orthogonal-to-Skin Method (POS)

Finally, we compared the RMSE values obtained from these methods.


### What is rPPG (Remote Photoplethysmography)?
rPPG is a method to estimate the heart rate (HR) of a person remotely without contact. Our project will estimate contactless HR using deep learning methods with a camera. This technology helps the patients to get rid of contact medical devices.


In this part, we used iPhys toolbox for implementations of CHROM, GREEN, ICA and POS methods.
#### iPhys: An Open  Non-Contact  Imaging-Based  Physiological  Measurement Toolbox
In the past few years a lot of attention has been given to methods for remotely measuring physiological signals using low-cost cameras.  Imaging PPG (iPPG) focuses on the measurement of volumetric changes in blood flow at distance from the body using imaging devices to capture changes in transmitted or reflected light. Imaging ballistocardiography (iBCG) typically leverages optical flow estimation to track the vertical motion of the head or body from a video sequence. Both iPPG and iBCG methods can be used to recover human vital signals.

This toolbox contains MATLAB implementations of a number of algorithms for non-contact physiological measurement. This will enable researchers to present results on their datasets using standard public implementations of the baseline methods with all parameters known. The toolbox includes implementations of many of the most commonly used baseline methods for imaging photplethysmography (iPPG) and image ballistocardiography (iBCG).