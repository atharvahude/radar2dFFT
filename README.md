# radar2dFFT
This project is created as a part of Udacity Sensor Fusion Nano degree Program 

## 2D CFAR Process Implementation

### Overview

This README provides an explanation of the implementation steps for the 2D Constant False Alarm Rate (CFAR) process in MATLAB, including the selection of training and guard cells, the calculation of the threshold, and the suppression of non-thresholded cells at the edges.

### Implementation Steps for the 2D CFAR Process

1. **FMCW Waveform Generation**:
    - Define radar specifications: frequency of operation, max range, range resolution, and max velocity.
    - Design the FMCW waveform with parameters like sweep bandwidth and chirp time.
    - Generate the transmit and receive signals for the moving target.

2. **Range Measurement**:
    - Reshape the beat signal into a 2D matrix.
    - Perform a 1D FFT on the beat signal to obtain the range information.

3. **Range Doppler Response**:
    - Perform a 2D FFT on the reshaped beat signal to generate a range-Doppler map.
    - Visualize the range-Doppler map using the `surf` function.

4. **CFAR Implementation**:
    - Select the number of training and guard cells in both dimensions.
    - Define an offset value for the threshold.

### Selection of Training, Guard Cells, and Offset

### To ensure the final result shows a distinct and clear peak, the guard cells are large enough to exclude the peak from the noise calculations.

1. **Training Cells**:
    - Training cells are used to estimate the noise level in the vicinity of the Cell Under Test (CUT). In this implementation, 8 training cells are selected in both the range and Doppler dimensions (`Tr = 8`, `Td = 8`).

2. **Guard Cells**:
    - Guard cells surround the CUT to prevent the signal from the CUT from influencing the noise estimate. In this implementation, 3 guard cells are selected in both dimensions (`Gr = 3`, `Gd = 3`).

3. **Offset**:
    - The offset is added to the average noise level to set a threshold for detection. In this implementation, an offset of 5 dB is used (`offset = 5`).

### Suppressing Non-thresholded Cells at the Edges

To suppress the non-thresholded cells at the edges of the range-Doppler map:

1. **Edge Handling**:
    - The CFAR process does not threshold cells at the edges due to the lack of sufficient training and guard cells. These edge cells are automatically set to zero in the thresholded output.

2. **CFAR Processing Loop**:
    - The loop iterates over the range-Doppler map, excluding the edges where training and guard cells cannot be applied.
    - For each CUT, the noise level is calculated from the surrounding training cells, excluding the guard cells.
    - The noise threshold is determined by converting the summed noise level to linear, averaging, converting back to logarithmic, and adding the offset.
    - The CUT value is compared to the noise threshold to decide if it should be marked as a detection (`1`) or not (`0`).