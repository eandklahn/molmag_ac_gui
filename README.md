# Molmag AC GUI

Import `process_ac` to use in scripts with your own data.

Install with
```
python -m pip install molmag_ac_gui
```
Run with
```
python -m molmag_ac_gui
```

To use the GUI for fitting relaxation times, your own file with data of (T, tau) can be loaded and fitted in the "Analysis"-tab
when the file is formatted as (the header line is mandatory, but the content is of no importance)

```
Temp;Tau
T1;tau1(;dtau1)
T2;tau2(;dtau2)
T3;tau3(;dtau3)

...
```
