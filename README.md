
<p align="center">
    <img src="PROSPECT.png" width="200">
</p>

---

<h1 align="center">
    <br>
    Positron Mobility Analysis
    <br>
</h1>

<h2 align="center">
    Calculating the displacement of the positron from the annihilation point. Expected result is ~0.05 cm away from the neutrino direction (Reactor). <a href="https://arxiv.org/pdf/hep-ex/9906011.pdf">Source:</a> Page 3, Section 2.1: Positron displacement.
</h2>

<h2>
    Summary
</h2>

The analysis code first creates two dataframes to read the .h5 file for the PG4 run. It inputs positron annihilation data (currently labeled NCapt), and the data of the primaries (Prim). All events that are not in positron annihilation are dropped from the primaries, as well as any events that are abnormal when it comes to capture time, annihilation products and energies, as well as segment numbers. These two dataframes are then combined and the difference between the *x*, *y*, and *z* values are calculated and averaged into a total vector. The two tables consider different axes as *x*, *y*, and *z* as well as the boundaries so cuts and offsets are applied. 

<h2>
    Requirements
</h2>

* Python 3+
    * pandas
    * numpy
    * matplotlib
    * h5py

<h2>
    How To Use (Debian)
</h2>

```bash
# Clone this repository
git clone https://github.com/shayshunk/PositronMobilityAnalysis

# Go into the repository
cd PositronMobilityAnalysis

# Install dependencies
sudo apt install python3-pandas python3-matplotlib python3-numpy python3-h5py

# Run the code. Make sure you have the .h5 file in this directory and edit line 6 to your filename.
python3 PostironAnnihilation.py
```

---
