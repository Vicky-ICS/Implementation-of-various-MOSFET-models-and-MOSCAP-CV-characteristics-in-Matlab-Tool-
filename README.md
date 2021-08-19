A code to generate HFCV and LFCV curves for NMOS and PMOS capacitors. Input variables: Gate workfunction, oxide thickness, channel doping and type (NMOS or PMOS). 

Generated CV curves for different MOS parameters .

Another code to extract MOS capacitor parameters from the CV curve that I generate from the above code. 

Modified the earlier developed CV code (HFCV and LFCV) to include effects of fixed charges (positive and negative, do separately) and interface traps. Consider fixed oxide charge density of 2x1012 /cm2 and interface trap density of 3x1012 /cm2 (both located at the Si/SiO2 interface). Use TOX = 2.2 nm and NA = 3*1017 / cm3. 

Numerically solved the Pao-Sah and Brews models for the following NMOSFET:  L = 1mm, TOX = 10nm, n+ poly-Si gate, selected NA for target VT of 0.8V, range of VG and VD sweeps is 0 – 5V.  

Ploted transfer (ID – VG) characteristics for different VD values and output (ID – VD) characteristics for different VG > VT values. Compared results from Pao-Sah and Brews models, and also with piece-wise equations (above threshold and sub threshold). Compared above threshold and sub threshold full and approximate expressions.

Re-did the above for a PMOSFET, with p+ poly-Si gate.

Used constant effective mobility values of 200 cm2/Vs (electrons) and 100 cm2/Vs (holes) for my calculations

