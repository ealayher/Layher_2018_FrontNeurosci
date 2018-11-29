# Layher_2018_FrontNeurosci
Processing code for the Frontiers in Neuroscience research article [Failure to Affect Decision Criteria During Recognition Memory With Continuous Theta Burst Stimulation](https://doi.org/10.3389/fnins.2018.00705): [**OSF**](https://osf.io/r73xg/)  

**DATA FILES**:   
[**tms1_sdt.csv**](https://github.com/UCSBMemoryLab/Layher_2018_FrontNeurosci/blob/master/tms1_sdt.csv)    
Normalized c and d' values from Study 1

[**tms2_sdt.csv**](https://github.com/UCSBMemoryLab/Layher_2018_FrontNeurosci/blob/master/tms2_sdt.csv)   
Normalized c and d' values from Study 2

[**tms_behave.xlsx**](https://github.com/UCSBMemoryLab/Layher_2018_FrontNeurosci/blob/master/tms_behave.xlsx)    
Spreadsheet of responses and signal detection values from Study 1 and 2

**CODE FILES**:   
[**normalized_criterion.m**](https://github.com/UCSBMemoryLab/signal_detection_theory/blob/master/normalized_criterion.m)   
Function that calculates the normalized criterion by residualizing <i>c</i> against <i>d'</i>.

[**sdt_measures.m**](https://github.com/UCSBMemoryLab/signal_detection_theory/blob/master/sdt_measures.m)   
Function that calculates discriminability (<i>d'</i>) and criterion values (<i>c</i>).
 
[**lmm_tms.m**](https://github.com/UCSBMemoryLab/Layher_2018_FrontNeurosci/blob/master/lmm_tms.R)   
Analyses data using linear mixed models and outputs figures.

[**summarizeWithin.m**](https://github.com/UCSBMemoryLab/Layher_2018_FrontNeurosci/blob/master/summarizeWithin.R)  
Function used by [**lmm_tms.m**](https://github.com/UCSBMemoryLab/signal_detection_theory/blob/master/lmm_tms.m) to calculate mean and standard deviation.
