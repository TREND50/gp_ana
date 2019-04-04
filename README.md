# gp_ana
Analysis tool for GP data (Rust)

Relevant tools are available with 3 script at present:

## [ReadData.py](https://github.com/TREND50/gp_ana/blob/master/readData.py)
Script to read GP35 data produced by the RUST DAQ software and manipulate it. Relies on ```pyef``` package to access data. Main functions are:

### display_events():
Loads data and displays timetraces 

### get_time():
Loads data, build trig time info and orders it in increasing order.

### def build_coincs():
Looks for causal coincidences between antennas (defined as Delta_trig < Delta_Pos/c0) and save these coincident events to file RRunId_coinctable.txt.
