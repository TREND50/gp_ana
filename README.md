# gp_ana
Analysis tool for GP data (Rust)

Relevant tools are available with 3 script at present:

## [readData.py](https://github.com/TREND50/gp_ana/blob/master/readData.py)
Script to read GP35 data produced by the RUST DAQ software and manipulate it. Relies on ```pyef``` package to access data. Main functions are:

### display_events():
Loads data and displays timetraces 

### get_time():
Loads data, build trig time info and orders it in increasing order.

### def build_coincs():
Looks for causal coincidences between antennas (defined as Delta_trig < Delta_Pos/c0) and save these coincident events to file RRunId_coinctable.txt.

## [reaRecons.py](https://github.com/TREND50/gp_ana/blob/master/readRecons.py)
Script to plot the results of spherical and plane reconstructions using files RRunId_planerecons.txt & RRunId_sphrecons.txt produced with the [gp_recons](https://github.com/TREND50/gp_recons) software. Main functions are:

### plot_delays(runid, coincid,...)
Computes antenna trigger time expected from reconstructed wave and plots these vs measured ones fro coincidence coincid in run runid. Computes associated Chi2. Perfect reconstruction should correspond to distribution along 1st bissector and Chi2 = 0. Read [TREND 2011 paper](https://arxiv.org/abs/1007.4359) for more details (see Fig 6 in particular).

### plot_recons(runid)
Plots various distributions of reconstructed events in run runid.

