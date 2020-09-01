# Sound Speed Profile Calculator
This repository provides tools for calculating the Sound Speed Profile from OOI Data stored in CSVs

## Running Script
To get the sound speed profile for the Oregon Slope Base OOI location for July 14-15 2019, start python in your terminal and run the following:

```python
import CTD_Tools.py
ssp = CTD_Tools.get_ssp_PC01A()
```