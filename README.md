# gp_ana
Analysis tool for GP35 data (recorded with Rust DAQ)

Relevant tools are available within 3 scripts at present:

## [readData.py](https://github.com/TREND50/gp_ana/blob/master/readData.py)
Script to read GP35 data produced by the RUST DAQ software and manipulate it. Relies on [```pyef```](https://github.com/TREND50/pyef) package to access data. Main functions are:

### display_events():
Loads data and displays timetraces. Note that in the present version the value plotted is 10^(signal) to correct for teh logarithmic amplification of the Power Detector present in the analog chain of the Electronic Board (see [the GRAND manual](https://github.com/TREND50/GRANDproto_manual/blob/master/manual.pdf) for details).   

### get_time():
Loads data, build trig time info and orders it in increasing order.

### build_coincs():
Looks for causal coincidences between antennas (defined as Delta_trig < Delta_Pos/c0) and save these coincident events to file RRunId_coinctable.txt, to be used for source reconstruction (see [```gp_recons```](https://github.com/TREND50/gp_recons) for details).

## [readRecons.py](https://github.com/TREND50/gp_ana/blob/master/readRecons.py)
Script to plot the results of spherical and plane reconstructions using files RRunId_planerecons.txt & RRunId_sphrecons.txt produced with the [```gp_recons```](https://github.com/TREND50/gp_recons) software. Main functions are:

### plot_delays(runid, coincid,...)
Computes antenna trigger times expected from reconstructed wave and plots these vs measured ones for coincidence coincid in run runid. Computes associated Chi2. Perfect reconstruction should correspond to distribution along 1st bissector and Chi2 = 0. Read [TREND 2011 paper](https://arxiv.org/abs/1007.4359) for more details (see Fig. 6 in particular).

### plot_recons(runid)
Plots various distributions of reconstructed events in run runid.

## [anaSLC.py](https://github.com/TREND50/gp_ana/blob/master/anaSLC.py)
Script to analyze GP35 slow control data saved in SRunID.yaml files. Main functions are:

### loopSLCEvents(boardID,RUNID):
Extracts SLC info from unit boardID in file SRunID.yaml, reduce them and saves result into file SLC_bBoardID.txt if they are more recent than latest data stored in that file.

### displaySLC(boardID,RUNID):
Displays content of SLC_bBoardID.txt results file within a time window hardcoded (berk) in the script.

# ToDo
- Develop statistical analysis of time pulses: mean of baseline, std of baseline, amplitude distribution, distribution of trigger time positions are interesting variables. This could be based on loopEvents() function in (obsolete) (script [GRANDproto_python/anaData.py](https://github.com/TREND50/GRANDproto_python/blob/master/anaData.py). Additionnaly, the rate of transient signals with random triggers, as well as the corresponding trigger positions, should be studied. See item 4.1.6 in [the GRAND manual](https://github.com/TREND50/GRANDproto_manual/blob/master/manual.pdf) for details on random data.
- Optimise coincidence search by using "true" distance between detection units (i,j) instead of max(Distance) as done presently (see line #118 in [readData.py](https://github.com/TREND50/gp_ana/blob/master/readData.py)).
- Study distribution of Chi2 as a function of unit ID in [readRecons.py](https://github.com/TREND50/gp_ana/blob/master/readRecons.py) in order to identify possible offsets in time tag and/or antenna position. This will require fetching from RRunID_coinctable.txt the exact list of unit IDs participating in a given coincidence. 
- Set up a mechanism to check if trig time info is indeed present in each event and flag antennas without such info.
- Implement calibration analysis. Could be based on (obsolete) script [GRANDproto_python/anaCalib.py](https://github.com/TREND50/GRANDproto_python/blob/master/anaCalib.py).
- Test analysis for large data file (now only up to 1GB).
