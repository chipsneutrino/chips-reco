
# chips-reco

CHIPS reconstruction software. Contains code for displaying, analysing, reconstructing, and classifying events. A heavily modified version of [WCSimAnalysis](https://lbne.bnl.gov/svn/people/CommonReconstruction/WCSimAnalysis/)

## Building

You need an up-to-date version of ROOT6 and chips-sim installed. This is provided by the chips-env and chips-sim repositories.
The first time you run reco-setup.sh it will build chips-reco, any further changes will require a 'make' command

```
$ source reco-setup.sh
```

## Running the Event Display

```
$ source reco-setup.sh
$ evDisplay
```

## Cleaning Everything Up

To remove all artifacts and return to the base state run...

```
$ source reco-tidy.sh
```

## Input/Output

Input to code is WCSim ROOT classes, using WCSimRootGeom and WCSimRootEvent to build
instances of WCSimOutputTree.

## Running the code

Running of the code except for the event display is done via ROOT macros which load
in the compiled libraries.

For examples look in the "config/example" directory
This directory also includes a sample WCSim .root file for testing

```
$ root -l -x -q wc_trackfitter_el_mu.C  # electron + muon fit 
$ root -l -x -q wc_trackfitter_el.C  # electron fit 
$ root -l -x -q wc_trackfitter_mu.C  # muon fit 
$ root -l -x -q wc_trackfitter_piZero.C  # pi-zero fit 
$ root -l -x -q wc_trackfitter_cosmic.C  # cosmic fit 
$ root -l -x -q wc_train_read_PID.C  # run the PID
```

## Notes on the PID

There are 5 different types of event sets that need to be used within the PID...
- numu_ccqe
- numu_all (i.e. CCQE + CC non-QE + NC events in the correct proportions)
- nuel_ccqe
- nuel_all (i.e. CCQE + CC non-QE + NC events in the correct proportions)
- numu_nc events

The event sets used for training the two ANN's are...
elANNmu -> Signal: nuel_ccqe , Background: numu_ccqe
elANNnc -> Signal: nuel_ccqe , Background: numu_nc
We use ccqe events so that they are trained on the cleanest subset of events.

A preselection is applied to any events that go into the training or are used as samples for testing...
- nHits < 50
- recoE_el < 550
- recoE_mu < 550
- recoE_el > 4950
- recoE_mu > 4950
- vtxRho_el > 1100
- vtxZ_el > 50 from top
- vtxZ_el < 50 from bottom
- veto
- fracHitsDownstream < 0

i.e. some activity in the detector, the reco energy not railing, and a 1m fiducial volume cut for the electron fit,
and a cut to remove escaping events that the veto would identify, as well as addition cuts that act as checks that the events have run correctly.

You need ~20k events of each type to get the appropriate stats and ensure that the training converges.
A large proportion (~80%) of numu_nc events are removed by the cut, mainly due to their elFit and muFit reco energies being < 550Mev. 
This is a good thing, however, leads to some issues with having a lack of NC events to train the elANNnc on.