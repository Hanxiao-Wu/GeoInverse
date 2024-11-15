# Monte Carlo Joint Inversion Program

## Overview
This program, written in C++, is designed for geophysical inversion applications using a classical Monte Carlo inversion approach. It provides a framework to estimate subsurface structure parameters with multiple types of geophysical data, particularly useful for studying crustal and upper mantle structures.

## Features
- Support for Multiple Data Types:
	- Receiver Functions (RF): Can take stacked RF to fit its waveform or multiple RFs from different events to fit theirs arrivals (H-k stacking)
	- Rayleigh wave Dispersion Data: Include both phases and group velocity. Potentially can also take Love wave dispersion, but not tested yet.
	- Rayleigh wave H/V ratio (Ellipticity)
	- Rayleigh wave Local Amplification:
- Inverted Results:
	- Vs (fine 1-D structure)
	- Vp/Vs (can be layered structure instead of just a bulk value)
	- Density (currently not tested)
	- Moho Depth
	- Temperature (future feature)
	- Pressure (future feature)
> [!IMPORTANT]
> - While the GeoInverse program is designed to be highly flexible, allowing a wide range of parameters to participate in the inversion, the reliability of the results fundamentally depends on the **input data**, not the inversion methodology itself.
> - Users must carefully choose inversion parameters based on the sensitivity of their data. For instance, although the program permits inversion for detailed Vp/Vs structures, such results are only trustworthy if the input data have a depth-related sensitivity to Vp/Vs variations. Otherwise, the inversion might produce a "result," but that result is just 'garbage in, garbage out'. Always assess the sensitivity and quality of your input data before proceeding.

## Inversion Setup
To run the program, three main files are required:
1. **Inversion Control File**(`*.control`)
   - This file defines the overall configuration of the inversion process, such as the number of Monte Carlo searches, the iteration count per search, and the data source weights.
   - It acts as the main setup file that connects all components of the inversion.
2. **Model File**
   - This file describes detailed 1D model using a smaller number of parameter.
   - The model described in this file also serves as the center of the model space.
3. `in.para`**File**
   - This file, together with the model file, defines the model space. It specifies which parameter to perturb, its bounds (absolute or percentage), and the step size for Monte Carlo sampling.

### 1. Inversion Control File (`*.control`)
The `*.control` file is essential for configuring inversion parameters and Monte Carlo settings. It tells the code where to read the data from, the weights for each dataset, how many searches to perform, how many iterations to run for each search, and so on. Below is an [example control file](test.control) with explanations for each parameter:
```
model 2 tar.mod        # Number of basic layers; model file
para in.para           # Parameter file for setting up perturbed parameters

disp R 4 p template_p.dat g template_g.dat e template_e.dat a template_a.dat  
# Type of surface wave data:
#  - "R" for Rayleigh waves
#  - 4 indicates the type of surface wave data, with the following files for phase (p), group (g), ellipticity (e), and amplification (a)

rf 2.5 0.06 template_rf.dat  
# Gaussian parameter for receiver function data
# Ray parameter for receiver function data
# File name for receiver function data

rfweight 0.4           # Weight for receiver function data

hk tar_hk.lst 1 0      
# File listing paths to all individual receiver function files
# Number of discontinuities
# Index of the discontinuity, starting from 0

hkweight 0.3 0.4 0.3   # Weights for each phase (Ps, PpPs, PsPs+PpSs)

monol 0                # index of the group which if forced to be monotonicly increasing

Eweight 1.0 0.035      
# Weight for energy generated from each discontinuity
# Reference energy used for normalization

#model -1               # Number of models in each search; -1 for forward calculation
#search 30              # Number of Monte Carlo searches; -1 for prior sampling
                         # If #model = -1, this parameter is ignored

outdir tar tar          # Output directory and file name

end                     # End of control file
```

### 2.Model Setup File
The program is designed to derive a detailed 1D subsurface model (e.g., Vs, Vp, density). Since a model with hundreds of layers would require an impractical number of parameters for Monte Carlo inversion, we instead adopt a layered approach:
1. Basic Layers: The model is divided into several main sections ( refered as group, usually two or three) to represent primary divisions:
	- Two-Layer Example: Crust and uppermost mantle.
	- Three-Layer Example: Sediment, crystalline crust and uppermost mantle.
2. Finer Layers within Basic Layers: Each basic layer is further divided into a finer 1D model. This detailed layering is controlled by user-defined parameters in the model file. Each fine layer's properties are determined by interpolation (e.g., B-spline) using a few coefficients, typically no more than five. This approach ensures a smoother transition in model parameters and reduces the inversion complexity.
For further information, refer to [Shen et al., 2013](https://academic.oup.com/gji/article/192/2/807/580799), and [Wu et al., 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023JB027952).
#### Model File
The model file defines the configuration of the basic layers and their fine structure, supporting the generation of a detailed 1D model using a variety of interpolation approaches. Also the model described in the model file serves as the **reference model** whicn will be the center of the model space. Below is a breakdown of each column's purpose.
##### Column Description
1. Group Index (Column 1):
	- Represents the index of the "basic layer" or **group**, starting from 0.
2. Parameter Type (Column 2):
	- indicates the property type for this row:
		- `1`- Vs(km/s)
		- `2`- Vp/Vs
		- `3`- Density(g/cm^3)
		- `4`- Qs (Quality factor for Vs)
		- `5`- Qp (Quality factor for Vp)
		- `6`- Temperature
		- `7`- Pressure
		
3. Interpolation Approach (Column 3):
	- Specifies how to generate the detailed 1D model within the group:
		- `1`- Gradient interpolation (linear gradient from top to bottom of the group)
		- `2`- Layered model (using this way ,this group is essentially already a layered model instead of interpolated from some coefficients)
		- `3`- B-spline interpolation (smooth curve fitted to the parameters)
		- `4`- Bulk value for entire group (using this way, if no additional anomaly added, then the whole group will be just one layer)
		- `-1`- Water layer (using this way, then the Vs will be set to be 0)
		- `-2`- Ice layer (to be done)
		- `-3`- Brocher's empirical relationship, only applied to Vp or Density (using this way the Vp or density will be scaled from Vs using corresponding empirical relationship)
		- `-4`- Hacker's empirical relationship, only applied to mantle Density

4. Group Thickness (Column 4): 
	- Total thickness of the group in km
	- For Vs, Vp, and density, the group thickness should generally be the same since discontinuities in these properties are assumed to align within a group. If temperature and pressure are included in the future, groups may have differing thicknesses, but currently, group thickness should **remain consistent for all properties within the same group**.

5. Number of Parameters (Column 5):
	- Specifies the number of parameters provided to generate the detailed 1D model within the group.
	- For example, if this number is `n`, then the next `n` columns contain the corresponding parameter values.

6. Parameters for Interpolation (Column 6 to 6+n):
	- These are the actual parameter values used for interpolation within the group.
	- The parameters vary by interpolation approach:
		- For `Gradient`, **two** parameters are required: the top and bottom values of the property in this group (e.g., top and bottom Vs).
		- For `Layered`, the number of parameters should match the number of fine layers.
7. Number of Anomalies (column 7+n):
	- Optional. Specifies the number of additional anomalies to add to the smooth 1D model within the grou (Let's say the number is `m` )p.

8. Anomaly Parameters (column 8+n to 8+n+3m):
	- Each anomaly requires three parameters:
		- **Upper boundary depth**: described in percentage of group thickness
		- **Lower boundary depth**: described in percentage of group thickness
		- **Anomaly value**: only allowed to be a bulk value currently
	- This section will contain `3m` values (three per anomaly)

9. Number of Fine Layers:
	- The desired number of layers in the fine 1D model for this group. The thickness of each fine layer is calculated as the **group thickness divided by the number of layers**.

##### Example
You can find an example model file [here](test.mod).
The following is an explanation in a group-by-group format. Each group has 7 rows to describe its Vs, Vp(or VP/Vs), density, Qs, Qp, T, P, respectively.

***group 0*** Sedimentary layer
```
0 1 1 1.0 2 1.5 2.7 0 5 0.
0 2 -3 1.0 0 0 5
0 3 -3 1.0 0 0 5
0 4 -3 1.0 0 0 5
0 5 -3 1.0 0 0 5
0 6 4 1.0 1 300 0 5
0 7 4 1.0 1 400 0 5
```
- Line 1: `0 1 1 1.0 2 1.5 2.7 0 5 0.`
	- **Group Index (0)**: Identifies the sedimentary layer as the first group in this model.
	- **Property Type (1)**: This line describes the shear wave velocity (Vs).
	- **Interpolation Approach (1)**: Uses gradient interpolation from top to bottom within this group.
	- **Group Thickness (1.0)**: Specifies the thickness of the sedimentary layer as 1.0 kilometers.
	- **Number of Parameters (2)**: Two parameters are needed for gradient interpolation.
	- **Parameter (1.5 2.7)**: Top (1.5km/s) and bottom (2.7km/s) Vs values for linear gradient.
	- **Number of Anomalies (0)**: No additional anomalies are added within this layer.
	- **Number of Fine Layers (5)**: This group is  divided into 5 fine layers.
	- **Depth of the Top Boundary (0.)**: Only the first group needs this value.
- Lines 2–5: `0 2 -3 1.0 0 0 5` to `0 5 -3 1.0 0 0 5`
	- **Property Types (2, 3, 4, 5)**: Define Vp/Vs, density, Qs, and Qp respectively.
	- **Interpolation Approach (-3)**: Uses Brocher’s empirical formula for automatic calculation based on Vs.
	- **Parameters**: No additional parameters are needed.
	- **Fine Layers (5)**: Matches the fine layer count of Line 1.
- Lines 6–7: `0 6 4 1.0 1 300 0 5` and `0 7 4 1.0 1 400 0 5`
	- **Property Types (6, 7)**: Define temperature (T) and pressure (P).
	- **Interpolation Approach (4)**: Uses bulk assignment with a single value across the layer.
	- **Parameter (300 for T, 400 for P)**: Sets temperature and pressure values.
	- **Fine Layers (5)**: Matches other properties in this group.
So basically, this group describes a 1.0 km thick sedimentary layer with a Vs gradient from 1.5 to 2.7 km/s, where the other properties of this layer (Vp, density, Qs, and Qp) are derived empirically from Vs. This way, these parameters can not be directly perturbed during the inversion. In the current version, the algorith doesn't do anything with the T and P properties, so 

**group 1** Crystalline Crust
- Line 8: `1 1 3 29.0 5 3.0 3.2 3.5 3.7 3.9 0 25`
	- **Group Index (1)**: Identifies the crystalline crust layer as the second group in this model.
	- **Property Type (1)**: This line defines Vs.
	- **Interpolation Approach (3)**: Uses B-spline interpolation for a smooth Vs profile.
	- **Group Thickness (29.0)**: Specifies thickness as 29.0 km.
	- **Number of Parameters (5)**: Five parameters for the B-spline.
	- **Parameters (3.0, 3.2, 3.5, 3.7, 3.9)**: B-spline control points for Vs.
	- **Fine Layers (25)**: Divides this group into 25 fine layers.

- Line 9: `1 2 4 29.0 1 1.71 1 0.5 1.0 1.78 25`
	- **Property Type (2)**: Represents Vp/Vs.
	- **Interpolation Approach (4)**: Bulk Vp/Vs across the group.
	- **Parameter (1.71)**: Vp/Vs ratio.
	- **Number of Anomalies (1)**: One anomaly is added.
	- **Anomaly (0.5 1.0 1.78)**:
		- **Upper Boundary (0.5)**: 50% depth of this layer.
		- **Lower Boundary (1.0)**: 100% depth of this layer.
	    - **Anomaly Value (1.78)**: Vp/Vs within anomaly.
	- **Fine Layers (25)**: Consistent with Vs layering.

- Lines 10–14: `1 3 -3 29.0 0 0 25` to `1 7 4 29.0 1 1000 0 25`
	- Define density, Qs, Qp (using empirical relations), temperature (700), and pressure (1000) in the crystalline crust.

So, layer 1 represents the 29 km thick crystalline crust, where Vs is smoothly interpolated with B-spline control points and the Vp/Vs structure is essentially divided into upper crust (top 50% of the crust) and lower crust (bottom 50% of the crust). Temperature and pressure are constant across the layer, while other properties are computed based on Vs.

**group 2** Uppermost Mantle
- Line 15: `2 1 3 150.0 5 4.2 4.35 4.45 4.53 4.6 0 30`
	- **Group Index (2)**: Identifies theuppermost mantle layer as the third group in this model.
	- **Property Type (1)**: Defines Vs.
	- **Interpolation Approach (3)**: B-spline interpolation for Vs.
	- **Group Thickness (150.0)**: Thickness of 150 km.
	- **Parameters (4.2, 4.35, 4.45, 4.53, 4.6)**: B-spline control points for Vs.
	- **Fine Layers (30)**: Divides this group into 30 fine layers.

- Line 16: `2 2 4 150.0 1 1.789 0 30`
	- **Property Type (2)**: Vp/Vs with bulk application.
	- **Parameter (1.789)**: Sets Vp/Vs ratio.

- Lines 17–21: `2 3 -3 150.0 0 0 30` to `2 7 4 150.0 1 1000 0 30`
	- Define density, Qs, Qp (empirical relations), temperature (800), and pressure (1000) in the mantle.

**Overall Summary**
The model file defines a three-layer structure (sedimentary layer, crystalline crust, and uppermost mantle), using a combination of linear, bulk, and B-spline interpolations. Each layer is characterized by specific Vs, Vp/Vs, density, Qs, and Qp, enabling a detailed 1D model suitable for forward calculation. 

> **Note:** The designations "sedimentary layer," "crystalline crust layer," and "uppermost mantle layer" are conceptual labels applied for user interpretation. The GeoInverse program does not explicitly recognize these as specific geological layers; rather, it organizes model groups from shallow to deep solely based on their `group index`. This indexing enables flexibility in setting up layers without requiring geological definitions in the input files.

### 3. `in.para` File
The `in.para` file is used to define the parameters for the Monte Carlo inversion, specifying which parameters will be perturbed and their perturbation range. Each row in this file corresponds to a single parameter for a specific group and describes how this parameter will be handled during the inversion process. The columns in this file serve different purposes, from identifying the group and parameter type to specifying boundary conditions and anomaly values.

#### Column descriptions
1. **Group Index** (Column 1): This column  which group the parameter belongs to. The index is based on the model's group numbering.
2. **Property Type** (Column 2): This column defines the property being described:
	- `0`: Thickness (km)
	- `1`: Vs (km/s)
	- `2`: Vp/Vs
	- ... (Basically the same as the model setting file)
 	- `-10`,`-11`,`-12`: Additional anomaly for Vs (top boundary, bottom boundary, value, respectively)
  	-  `-20`,`-21`,`-22`: Additional anomaly for Vp/Vs (top boundary, bottom boundary, value, respectively)
   	- ... (Similarly patterns for other anomalies) 
3. **Absolute/Percentage Value Indicator** (Column 3): This column defines whether the value in the fourth column is an absolute value or a percentage of the model space.
	- `0`: The value in the fourth column is expressed as a **percentage** of a reference value (which is the value in the model file).
	- `1`: The value in the fourth column is an **absolute** value. 
4. **Model Space Radius** (Column 4): This column defines the radius (or model space variation) of the parameter in the specific dimension. This represents the extent to which the parameter can vary during the inversion process.
   - The final model space will be the value in the model file plus/minus this radius
5. **Monte Carlo Step Size** (Column 5): This column defines the step size for the Monte Carlo inversion process, representing the standard deviation of the Gaussian distribution from which perturbations are drawn. This determines how much the parameter can change in each iteration of the Monte Carlo search.
6. **Parameter index** (Column 6): This column is used to identify the specific parameter when there are multiple parameters of the same type within the same group.
   - The value indicates the **sequence number** of the parameter within the group. For example, if there are two Vs values in group 0, the first one would be assigned `0` and the second one would be assigned `1` in this column.

#### Example
[Here](in.para) is an example `in.para` file, with each row explained:
- Line 1: `0 0 0 1.0 0.1`
	- This row defines the thickness (telled by the 2nd number `0` in this row) of group 0 (telled by the 1st number `0`)
 	- The radius of this model space is 100% (the 3rd number `0` tells us the radius is defined in percentage, the 3th number `1.0` means 100%) of the reference value (which is 1.km - set in the model file)
  	- It has a step size of `0.1` (telled by the 5th number) for MC interations.
  	- Thickness is described by a single value in one group so it does not need the 6th column.
  	- So in summary, the thickness of the group 0 is perturbed between 0.0~2km, with a step size of 0.1 km.
- Line 2: `0 1 1 0.5 0.05 0`
	- This row defines the first (telled by the 6th number `0`) Vs (telled by the 2nd number `1`) parameter for group 0 (telled by the 1st number `0`).
   	- The radius in model space is `0.5` (absolute value,telled by the 3rd and 4th number), with a step size of `0.05`.
- Line 3: `0 1 1 0.5 0.05 1`
  	- This row defines the second (telled by the 6th number `1`) Vs parameter for group 0.
  	- The rest is the same as Line2
- Line 4:`1 0 1 5.0 0.5`
  	- The thickness of the group 1 is perturbed between `29-5=24` to `29+5=34`km, with a step size of 0.5 km
- Line 5:`1 1 1 0.5 0.05 0`
  	- The **first** Vs parameter (which is a B-spline coefficient according to the model file) is perturbed between `3.0-0.5=2.5` and `3.0+0.5=3.5` km/s, with a step size of 0.05.
- Line 6: `1 1 1 0.5 0.05 1`
  	- The **second** Vs parameter (which is a B-spline coefficient according to the model file) is perturbed between `3.2-0.5=2.7` and `3.2+0.5=3.7` km/s, with a step size of 0.05.
- Line 7 to Line 9: ...(Similar to Line5, Line6)
- Line 10: `1 2 1 0.15 0.02 0`
  	- The first **Vp/Vs** parameter (telled by the 2nd number `2`) is perturbed between `1.71-0.15` and `1.71+0.15`, with a step size of 0.02
- Line 11: `1 -22 1 0.15 0.02 0`
  	- The first Vp/Vs **anomaly value** (telled by the 2nd number `-22`) is perturbed between `1.78-0.15` and `1.78+0.15` with a step size of 0.02
- Line 12 to Line 16: ... (Similar to Line5~9)

## Input Data Format
The **surface wave dispersion, H/V ratios, and waveform-fitting receiver functions** are stored in plain `.txt` files. The first row of the file specifies the number of rows and columns in the dataset (exclude the first row). From the second row onward, the data is structured as follows:
- Surface Wave Dispersion (km/s) and H/V Ratios:
	- Column 1: Period (s)
	- Column 2: Measured Value
	- Column 3: Measurement Error
- Surface Wave Local Amplification
	- Column 1: Period (s)
   	- Column 2: Measured Value (This is a relative value)
     	- Column 3: Measurement Error
	- Column 4: Reference Value
- Receiver Functions (For waveform fitting)
	- Column 1: Time (s)
   	- Column 2: Amplitude
     	- Column 3: Measurement Error
The **receiver functions that is used for H-kappa stacking** are stored in **SAC Format**.These SAC files include the following:
- Standard SAC headers (beggining time, number of data points, sampling rate, and so on).
- The user3 variable in the SAC header stores the ray parameter of the corresponding receiver function.
Some example data can be found [here](FowardTest) (The file whose name contains 'template')
