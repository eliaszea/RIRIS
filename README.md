RIRIS Room Impulse Response Interpolation with Shearlets
-----------------------------------------------------------------

**1. INTRODUCTION**

RIRIS (Room Impulse Response Interpolation with Shearlets) contains the MATLAB implementation of the algorithm in [^1], on RIR interpolation using shearlet dictionaries. 

[^1]: [E. Zea, “Compressed sensing of impulse responses in rooms of unknown properties and contents,” J. Sound Vib.  459, 114871 (2019)](http://kth.diva-portal.org/smash/record.jsf?pid=diva2%3A1340771&dswid=-7099).

RIRIS is covered by a GPL v3 license (see LICENSE for license terms).

**2. INSTALLATION**

Download the .zip file and extract it in your folder of preference. There is no need to install files, as the dependencies are added (and removed) automatically when running the main script. 

OBS! The folder called ‘dependencies’ must be in the same directory as the main script. The computer storage required by RIRIS is 170 MB. 

**2. MAIN SCRIPT**

_2.1. RIRIS_

Run this script to perform the iterative thresholding of shearlet coefficients in a room of choice (‘Munin’,’Freja’,’Balder’), with specific no. decomposition scales and no. iterations. 

The optimal threshold parameter is read from ‘dependencies/regularizationData’ so as to decrease the running times of the whole script. If the user wishes to compute the Pareto curves and estimate the optimal threshold anyway, then the binary variable ‘paretoFlag’ in Line97 should be set to 1. OBS: this may incur overall longer computation times of RIRIS. 

Once the interpolation is finished, Figure 1 outputs the results in the form of RIR images (under-sampled, interpolated, and reference), and Figure 2 outputs the modal assurance criterion (MAC) vs. frequency. In addition, the normalized mean-squared error (NMSE), as well as the computation time (CT), are displayed in the Command Window.

Example of usage: interpolate RIRs in lecture room 'Munin', provided an under-sampling ratio of 3, a shearlet dictionary with 4 scales, and running 20 thresholding (ISTA [3]) iterations:

```
RIRIS(‘Munin',3,4,20);
```

OBS! A fifth input argument (saveFlag) is accepted by RIRIS, which, if set to TRUE, stores the results (NMSE, MAC, image, image_recov,…) into a .mat file in the folder ‘dependencies/results.’

**3. DEPENDENCIES**

_3.1. BASIS FUNCTIONS_

Cone-adapted shearlets for every room (i.e., every $(T,M)$ combination) and for values of $\tau$ (decomposition scales) $2$ to $5$. 

Example: ‘Balder_tau_4.mat’ contains the $K \times T \times M$ shearlet bases used in the meeting room “Balder” provided a 4-scale dictionary (i.e. $K = 61$).

These basis functions are computed with the [Fast Finite Shearlet Transform toolbox](https://github.com/rujieyin/toolbox_FFST), copyright 2014 Sören Häuser.

_3.2. MEASUREMENT DATA_

Reference RIR images for each room. These _.mat_ files are loaded with the utility function ‘loadRIRs’. The content of each file is a structure ‘out’ with fields: 

	- image:  	T x M numeric array 
	- fs: 		temporal sampling frequency
	- T:		no. time samples
	- M: 		no. microphones
 
_3.3. REGULARIZATION DATA_

Optimal threshold parameters $\beta^\star$, for every room and no. decomposition scales, at which $J(\beta)$ attains a maximum (trade-off sparsity and model misfit). 

Example: _‘Freja_u3_tau3.mat’_ contains the value of $\beta^\star$ in Freja, given an under-sampling ratio of 3, and 3 decomposition scales. 

In the main script, the default is to load these parameters (avoiding the need to maximize the curvature every time). If one wants to run the maximization of $J(\beta)$ anyway, setting ‘paretoFlag’ to 1 in Line 89 does the trick. Be careful on where these end up stored.

The content of each .mat file is a structure ‘reguThresh’ with fields:
- beta_star: 	value of $\beta$ at which $J(\beta)$ attains a maximum
- Jcurve:	curvature function $J(\beta)$
- beta_set: 	pool of values of $\beta$ for maximization of $J(\beta)$

_3.4. UTILITY FUNCTIONS_

Below is an overview of the utility functions that are used in the main script. More detailed explanations can be found in their individual preambles. 

3.4.1. computePareto: compute L-curve corner with cubic spline interpolation

3.4.2. expansionMtrx: generate expansion matrix restricting vertical-like shearlets

3.4.3. extRIRimage: perform spatial extrapolation of RIR image

3.4.4. ffst: vectorized fast finite shearlet transform (both in/out are vectors)

3.4.5. iffst: same as 3.4.4. but inverse transform

3.4.6. ista: iterative soft thresholding algorithm

3.4.7. loadRIRs: load RIR reference image

3.4.8. lpbp1D: 1D linear predictive border padding of spatial dimensions

3.4.9. perforMetrics: assessment function that outputs performance metrics and figures

3.4.10. selectionMtrx: generate masking (selection) matrix


**4. RELEASE HISTORY**

	Release #1	 RIRIS v1.0 	E. Zea	2023-04-27


**5. FEEDBACK & CONTACT INFORMATION**

Your questions, suggestions, and feedback can help improve the quality of this software. Feel free to contact me at

	Elias Zea (zea@kth.se)
	Marcus Wallenberg Laboratory for Sound and Vibration Research
	KTH Royal Institute of Technology
	Teknikringen 8
	10044 Stockholm, SWEDEN


**6. LEGAL INFORMATION**

Copyright 2019 Elias Zea

This software was written by Elias Zea, and it was created during a postdoctoral period at the Marcus Wallenberg Laboratory for Sound and Vibration Research, KTH Royal Institute of Technology. 

RIRIS is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. If not stated otherwise, this applies to all files contained in this package and its sub-directories. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
