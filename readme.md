### Go to GitHub Pages:
 https://shri-amogh.github.io/Planetarium/

.
# ALL IMPORTANT MATH LISTED BELOW


## 1, time standardization

All astronomical calculations require converting civil time into uniform, measured from standard J2000 Epoch


count of days since noon on January 1, 4713 BC.

$$
JD = D + \lfloor \frac{153M' + 2}{5} \rfloor + 365Y' + \lfloor \frac{Y'}{4} \rfloor - \lfloor \frac{Y'}{100} \rfloor + \lfloor \frac{Y'}{400} \rfloor - 32045
$$



$T$ is the number of 36,525-day Julian centuries elapsed since the standard epoch J2000.0 (JD 2451545.0). This value is used for calculating time-dependent orbital elements.

$$
T = \frac{JD - 2451545.0}{36525}
$$

---

## 2, heliocentric coordinates 

The planet's position relative to the Sun

Each Keplerian element (e.g., Semi-major axis $a$, Eccentricity $e$, etc.) is modeled as a polynomial function of $T$:
$$
\text{Element} = \text{Element}_0 + \text{Element}_1 T + \text{Element}_2 T^2 + \dots
$$


#### mean anomaly

The angle representing the average position of the planet in its orbit. It is often calculated directly from the Mean Longitude ($L$), Longitude of Perihelion ($\varpi$), and Longitude of the Ascending Node ($\Omega$):

$$
M = L - \varpi \quad \text{where } \varpi = \omega + \Omega
$$


Due to the elliptical nature of the orbit, the true position is found using the Eccentric Anomaly ($E$), which is solved iteratively (e.g., using the Newton-Raphson method):

$$
M = E - e \cdot \sin(E)
$$


$E$ is converted to $v$ (the true angular position) and $r$ (the distance from the Sun):

$$
\tan\left(\frac{v}{2}\right) = \sqrt{\frac{1 + e}{1 - e}} \cdot \tan\left(\frac{E}{2}\right)
$$

$$
r = a (1 - e \cdot \cos(E))
$$

Heliocentric ecliptic rectangular coordinates ($x_{hel}, y_{hel}, z_{hel}$)
# Why rectangular coordinates? becasuse i dont like polar or spherical : )
 
The planet's position relative to the Sun, aligned with the Ecliptic plane. This often involves a series of rotations using the Argument of Perihelion ($\omega$), Inclination ($I$), and Longitude of the Ascending Node ($\Omega$):

$$\begin{aligned} x_{hel} &= r \cdot (\cos(\Omega) \cos(\omega + v) - \sin(\Omega) \sin(\omega + v) \cos(I)) \\ y_{hel} &= r \cdot (\sin(\Omega) \cos(\omega + v) + \cos(\Omega) \sin(\omega + v) \cos(I)) \\ z_{hel} &= r \cdot \sin(\omega + v) \sin(I) \end{aligned}$$


---

## 3, coordinate transform

align the coordinates with the Celestial Equator to find the observable position.



The Earth's position ($\vec{E}_{hel}$) is calculated using the same method as above. The planet's position relative to the Earth is the vector difference:

$$
\vec{P}_{geoc} = \vec{P}_{hel} - \vec{E}_{hel}
$$

$$\begin{aligned} X_{geoc} &= X_{hel} - X_{Earth} \\ Y_{geoc} &= Y_{hel} - Y_{Earth} \\ Z_{geoc} &= Z_{hel} - Z_{Earth} \end{aligned}$$



The coordinates are rotated from the Ecliptic plane to the Celestial Equator using the **Obliquity of the Ecliptic ($\epsilon$)**:

$$
\epsilon \approx 23.439291^\circ - 0.01300417^\circ T
$$

The rotation matrix $R_x(\epsilon)$ is applied:

$$
\begin{bmatrix} X_{eq} \\ Y_{eq} \\ Z_{eq} \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos(\epsilon) & -\sin(\epsilon) \\ 0 & \sin(\epsilon) & \cos(\epsilon) \end{bmatrix} \begin{bmatrix} X_{geoc} \\ Y_{geoc} \\ Z_{geoc} \end{bmatrix}
$$



The final output is the planet's angular position on the celestial sphere.

$$
R = \sqrt{X_{eq}^2 + Y_{eq}^2 + Z_{eq}^2} \quad \text{(Distance to Planet)}
$$

$$\begin{aligned} \delta &= \arcsin\left(\frac{Z_{eq}}{R}\right) \quad \text{(Declination)} \\ \alpha &= \arctan2(Y_{eq}, X_{eq}) \quad \text{(Right Ascension)} \end{aligned}$$
