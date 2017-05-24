# LactateThresholdAlgorithms
## Comparison of different algorithms/methods for lactate threshold detection.

### Methods:
* Baldari: second increase of at least 0.5 mmol/l
* Bishop: modified D-Max method
* Cheng: D-Max method
* Dennis: 3-segment linear regression
* Dickhuth: AT + 1.5 mmol/l
* Heck: 4 mmol/l-threshold
* Morton: 3-segment parabolic regression

### Example:

|         | Heck | Baldari | Dickhuth | SegLin | SegPara | Dmax | Dmod | mean | std |
|---------|------|---------|----------|--------|---------|------|------|------|-----|
| PB03    | 388  | 360     | 380      | 357    | 368     | 311  | 382  | 364  | 26  |
| PB05    | 283  | 260     | 272      | 253    | 236     | 225  | 276  | 258  | 21  |
| PB06    | 267  | 220     | 251      | 258    | 259     | 225  | 262  | 249  | 19  |
| PB08    | 215  | 180     | 203      | 219    | 214     | 194  | 224  | 207  | 16  |
| PB09    | 220  | 200     | 212      | 218    | 223     | 195  | 226  | 213  | 12  |
| Average | 275  | 244     | 263      | 261    | 260     | 230  | 274  | 258  | 16  |

![Example](https://raw.githubusercontent.com/A-A-G/LactateThresholdAlgorithms/master/lactate/LactateThreshold2_01.png)

### Required:
* MATLAB > 2016b (functions in script)
* Optimization Toolbox (fmincon)
* Bioinformatics Toolbox (suptitle, not necessary)

### Citation:
Send an [email](mailto:alexander@artigagonzalez.de) for more details.
