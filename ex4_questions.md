# Systems BioMedicine Exercise 4 - Lipidomics
#### Group 2

### 1. Lipidomics Analytics

##### a) Briefly describe the way of an analyte from injection to (fragment) detection in a Quadrupole - Orbitrap combination
The analyte is injected into the MS using electrospray ionization (ESI). It travels to the quadrupole, where molecules are trapped by switching electric fields at high frequency. Only ions with a specific m/z range can pass through. These are then fragmented using high- or low-pressure cells. The fragments are sent to the orbitrap, where they oscillate around the trap. Using a fourier-transformation, the fragment m/z can be calculated from the oscillation frequency. Using this second Spectrum, the exact molecule can be determined.

##### b) Shortly explain why hyphenation techniques like liquid chromatography help with the identification and quantification of lipids
If all molecules were injected into the MS at once, they would also be measured at the same time, making it difficult to distinguish molecules with a similar m/z ratio. Different molecules take a different amount of time to pass through a chromatograph, and are then measured in the MS at different times. Since we now have multiple spectra with time stamps, their peaks are no longer overlapping.

### 2. Spectra Processing

##### a) Briefly explain...

  - **...what a lockmass is and what it is used for**

    An empty sample is injected into the MS. Since there is nothing to analyze, the resulting spectrum shows the background caused e.g. by instrumental errors. This can then be used to correct for the background in the actual MS run.

  - **...how normalisation with internal standards work**

    By adding the same amount of internal standard to each sample, the peak intensity between samples, batches or experiments can be normalized, since the intensity of the internal standard should be consistent.

  - **...what peak centroiding is and how it works**

    A MS peak is defined through a gaussian fit applied over a number of data points. The centroid is the data point with the highest intensity, however it is usually not exactly equal to the highest point of the curve. The peak centroid is shifted, until it fits exactly onto the gaussian curve. Its' m/z and intensity values are therefore corrected.

##### b) What are isotopic patterns and why are they important to consider for feature identification

Something about isotopic peak clusters I guess?
