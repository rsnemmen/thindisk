Thin accretion disk spectrum
==============================

This python library computes the electromagnetic spectrum of a simple thin accretion disk, as a sum of local black body spectra. The units of the code are CGS. 

## How to use

Given the array `lognu` this method computes `log(nu*Lnu)` tabulated at the input frequencies for the given parameters.

```python
import thindisk as disk
import numpy as np

# array of frequencies
lognu=np.linspace(13,15,50)
# disk spectra
lum=disk.thindisk2(lognu,4e6,1e-3,3.,30.)
```

### Arguments to `thindisk2`

1. array of log10(frequencies)
2. black hole mass in solar masses
3. mass accretion rate in Eddington units
4. inner radius of the disk in Schwarzschild radii
5. inclination angle in degrees
6. outer radius in Schwarzschild radii, default is 1E5

The output is an array with `log10(nu*Lnu)`.


## More examples

See the accompanying jupyter notebook `spectra.ipynb`. Please inspect the source code for other useful methods.