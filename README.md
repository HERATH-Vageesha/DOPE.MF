# DOPE.MF
DNA Origami Protractor for Exploring Molecular Forces
TABLE OF CONTENT: 
01) Fiduciary markers and image Registration steps
02) Drift Correction for timelapses
03) Quantification of image intensity
04) Defining puncta uncertainy value
05) EMCCD imaging conditions
06)Defining End-to-End Distance, PHI and THETA Angles

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
