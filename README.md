
<p align="center">
    <img src="PROSPECT.png" width="200">
</p>

---

<h1 align="center">
    <br>
    (Anti)Neutrino Directionality Study
    <br>
</h1>

<h2 align="center">
    Calculating the average antineutrino direction in the PROSPECT detector. The source is the HFIR reactor. The reconstructed direction is achieved by tracking average delayed neutron displacement from the prompt positron in Cartesian before converting the coordinates to spherical. Paper in progress. 
</h2>

<h2>
    Summary
</h2>

The analysis code first reads analyzed PROSPECT data passed through the analysis framework PROSPECT2x and creates histograms to track the average displacement between the prompt and delayed events. The *x* and *y* positions are approximated as the midpoint of the segment in which the events occurred. The inoperative segments requires a more complex selection of events, which will be described in the paper. This also leads to a more complex error calculation which is outlined in the code. After the events are counted, the background subtraction is carried to get the IBD signal events which are then put through the error calculations and angle reconstructions. 

<h2>
    Requirements
</h2>

* ROOT 5+
* C++ compliler (I'm using g++ in the instructions)
* Analyzed and calibrated PROSPECT data 

<h2>
    How To Use (Debian)
</h2>

```bash
# Clone this repository
git clone https://github.com/shayshunk/NeutrinoDirectionality

# Go into the repository
cd NeutrinoDirectionality

# Install dependencies
sudo snap install root-framework

# Run the code. Make sure you edit lines 199 and 219 to point to your PROSPECT data files
g++ DeadSegmentCorrectionv5.cc -o DeadSeg `root-config --cflags --glibs`
./DeadSeg
g++ NeutrinoDirectionalityPlots.cc -o Plots `root-config --cflags --glibs`
./Plots
```

---
