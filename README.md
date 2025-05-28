# NTKSPH
NTK's 3D Speaker Scanning Spherical Harmonics Toolbox

## Introduction
This code is heavily inspired on the excellent work by NTK over on [Audio Science Review](https://www.audiosciencereview.com/forum/index.php?threads/diy-3d-speaker-scanner-the-mathematics-and-everything-else.9970/). It allows Sound Field Separation as measured by a Near Field Scanner. Look here for more info:

- [Near Field Scanner code](https://github.com/TomKamphuys/NFS)
- [Build thread](https://www.diyaudio.com/community/threads/klippel-near-field-scanner-on-a-shoestring.318151/)

The code is developed using [Octave](https://www.octave.org/), a [Matlab](https://www.mathworks.com/products/matlab.html) clone. It should also run in Matlab or do so with few adjustments, but I've never actually tried it.

## Installation
NTKSPH uses a few other projects to aid in the steps to process the data.

- [ita toolbox](http://www.ita-toolbox.org) (don't know if it is actually used). Read its README.md for how to find more info
- [MATAA (Mat's Audio Analyzer)](https://github.com/mbrennwa/mataa)
- some of maiky76's code he published somewhere diyadio or audiosciencerevies (don't know if it is actively used at the moment)
- ...


## Setup

Add NTKSPH and it subdirectories to the path.
Add mataa and its subdirectories to the path.


## Commands

Make sure you are in the NTKSPH directory.

### Process the data into a handy format for Octave/Matlab

```
[p, r, theta, z, f] = read_nfs_measurements('21062024-spherical');
save -7 21062024_measurement.mat r theta z p f
```

21062024-spherical is the directory with all the measurements (.wav files) as taken by the Near Field Scanner. This directory has to be present
in the Measurements subdirectory of NTKSPH.

### Do the actual fitting and Sound Field Separation

This step does all the heavy work on might take some time. A waitbar will be shown so you have something to look at :-) 

```
[CD, fit_error, CD_tot] = convert_to_coefficients();
```

CD contains the coefficients describing the outgoing sound
fit_error is the ... fit error
CD_tot contains the coefficients describing the in- and outgoing sound

#### Have a look at the fit error

```
figure
plot(fit_error)
semilogx(f, fit_error)
semilogx(f, fit_error)
xlabel('Frequency [Hz]')
ylabel('Fit error [dB]')
```

### Save data for quick reuse later

```
save -binary 22082024_CD8_acoustic_center.mat CD
```

22082024_CD8_acoustic_center.mat is the filename
CD is the variable to save (the outgoing sound coefficients in this example)

### Make 2D directivity plot

Command for vertical directivity:

```
[out, angle, freqs] = vertical_directivity(CD, f);
```

Command for horizontal directivity:

```
[out, angle, freqs] = horizontal_directivity(CD, f);
```

Now you can plot the data:

```
figure
pcolor(freqs, angle, out)
shading flat
set(gca,'xscale','log');
colormap jet
colorbar
max(out(:))
caxis([-70 -30]) % adjust for actual limits
```

### 2034CEA (Spinorama) Plot

TODO