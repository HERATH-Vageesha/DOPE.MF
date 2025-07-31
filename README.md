# DOPE.MF
DNA Origami Protractor for Exploring Molecular Forces
TABLE OF CONTENT: 
01) Fiduciary markers and image Registration steps
02) Drift Correction for timelapses
03) Quantification of image intensity
04) Defining puncta uncertainy value
05) EMCCD imaging conditions
06) Defining the PHI and THETA Angle



CONTENT

01. Fiduciary markers and image Registration steps


Fiducial markers which are diffraction limited with high photostability (TetraSpec Beads TM) were chosen for the procress of image registration. 
Further product information can be found below:
https://www.thermofisher.com/order/catalog/product/T7279

Imagine registration was performed via two seperate plugins available on FIJi/ImageJ namely



a) TurboReg

       Open ImageJ/Fiji --> Plugin --> TurboReg --> Affine Registration/Accurate/Automatic


b) Registration Plugin

      ImageJ/Fiji --> PLugin --> Registration --> Descriptor Based Registration (2d/3d)



Once the image has been registered, the image registration matrix is calulted. Using this matrix the experimental timelapses can be automatically registered moving forward.

02) Drift Correction for timelapses

To perform drift correct, both Picasso and Fiji/ImagJ can be used

ImageJ/Fiji: 

    ImageJ/Fiji --> Plugin --> StackReg --> Affine

Picasso: 

Please reffer to:
https://picassosr.readthedocs.io/en/latest/


03) Quantification of image intensity


Image intensity quantification is done via ThunderStorm Analysis which is a plugin available through Fiji/ImageJ. For more information please refer to

M. Ovesný, P. Křížek, J. Borkovec, Z. Švindrych, G. M. Hagen. ThunderSTORM: a comprehensive ImageJ plugin for PALM and STORM data analysis and super-resolution imaging. Bioinformatics 30(16):2389-2390, 2014.


04) Defining puncta uncertainy value

For more information please visit

Mortensen, Kim I., et al. "Optimized localization analysis for single-molecule tracking and super-resolution microscopy." Nature methods 7.5 (2010): 377-381.
Given then end to end distance of the probe, puncta with an uncertainty value lesser than 5nm were chosen and filtered when moving onto identifying probe angles. Given the uncertainty of either dipole being able to have a combine 10nm uncertainy would be the maximum uncertainy.

05) EMCCD imaging conditions


Pixel size: 150 nm


Photoelectorns per A/D count: 15.1


Quantum Efficiency; 1.0


Base levele [A/D counts: 166


EM Gain: 800


Peadout Noise [e-/pixel}: N/A


7) Defining the PHI and THETA Angle

The PHI angle was the orientation/dipole direction from the cy3b channel (TIRF 560) to ATTI647 channe; (TIRF 647) channel on a 2D plane when visualizing the angles from above. 
