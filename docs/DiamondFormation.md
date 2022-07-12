# Diamond formation
The diamond formation is estimated from the XRD signal by fitting peaks to the spectrum.


<span style="color:blue">TODO</span>
```
All functions for these calculations should be put to individual files and the notebooks only be used for the plotting! 
```

## Diamond density from peak position
Bragg equation leads the lattice plane distance $d$:
$$
n\lambda = 2d\sin(\theta) \\
d = \frac{n\lambda}{2\sin(\theta)}
$$
which is related to the lattice paramter $a$ via the Miller indices $h$, $k$ and $l$ via:

$$
a = d \cdot \sqrt{h^2 + k^2 + l^2} 
$$

therefore the volume of a unit cell is given by:
$$
V = a^3 = \left[\sqrt{h^2 + k^2 + l^2}\frac{n\lambda}{2\sin(\theta)}\right]^3
$$

The mass of a unit cell $m_{uc}$ is given by the number of atoms per unit cell $n_{uc}$ and the molar mass of species $x$ $m_x$.

$$
m_{uc} [kg] = n_{uc}(m_{x}*u)
$$

with $u = 1.66053904\times10^{-27}$. Multiplied by $1000$ to get the mass in $g$:

$$
m_{uc} [g] = 1000 * n_{uc}(m_{x}*u)
$$

Therefore the density is given by: 

$$
\rho = \frac{m}{V} = \frac{1000 * n_{uc}(m_{x}*u)}{\left[\sqrt{h^2 + k^2 + l^2}\frac{n\lambda}{2\sin(\theta)}\right]^3}
$$

the uncertainty can be calculated by: 

$$
\sigma_{\rho}^2 = \left|\frac{\partial\rho}{\partial\theta}\right|^2 \sigma_\theta^2 = \left| \frac{1000 * n_{uc}(m_{x}*u)}{\left[\sqrt{h^2 + k^2 + l^2}\right]^3}\frac{8}{\left(n\lambda\right)^3}\cdot3\cdot\sin^2(\theta)\cos(\theta)\right|^2 \cdot \sigma_\theta^2
$$

## Crystalite size from modified Scherrer equation given [here](https://cyberleninka.org/article/n/1113245.pdf)

In the Scherrer equation the crystalite size is related to the peak broadening by:
$$
L = \frac{K\lambda}{\beta \cos(\theta)}
$$
where $\beta$ is the full width at half maximim (FWHM) of the peak.
To combine the information from multiple peaks, the $\ln(\beta)$ is plotted against $\ln(\frac{1}{\cos(\theta)})$. By fitting a linear funtion trough the given points the crystalite size is given by the interception with the $y$ axis as $f(x=0) = \ln(\frac{K\lambda}{L})$. Hereby the $K$ is a factor which changes with crystalite shape and $\lambda$ is the X-ray wavelength. 