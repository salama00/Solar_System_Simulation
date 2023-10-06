# What the Project is About

This project is based on the paper "Migration and dynamical relaxation in crowded systems of giant
planets" by Adams & Laughlin. It aims to reproduce the results of the paper.
In a talk I discussed the main part of the paper: simulations of systems with 10 giant planets ranging
from 5AU to 30AU in circular uninclined orbits for 4 different mass distributions. The mass distributions
are: all planets have mass mj , all planets have mass 2mj , planets have masses from a uniform distribution of
the range 0 < mp < 4mj , and of the range −1 < log(mp/mj) < 1 (where mj is Jupiter’s mass and mp is the
planet’s mass). In addition, the author did a simulation in which a system contained one large planet and
19 smaller bodies. The large planet has mass equal to Jupiter’s and has equal semi-major axis (5AU). The
smaller bodies were randomly placed in the range 3 AU < r < 30 AU with masses randomly distributed in
the range 0 < mp < 0.5mj . Similar to the previous ensembles, 100 simulations were run. Seventy percent of
the time, the large planet becomes the innermost planet. When the large planet is the innermost planet its
orbital elements (on average) are ⟨a⟩ = 2.08AU ± 1.28AU and ⟨ϵ⟩ = 0.36 ± 0.20 (Adams & Laughlin 2003).
The average eccentricity is smaller compared to the scenario where there are 10 giant planets. This is most
likely because the scattering interactions due to the small bodies are likewise smaller (recall that a scattering
interaction takes away energy E = −GMmp/2r which is proportional to the mass). Despite the implication
of the uncertainty of the semi-major axis (a = [0.8AU, 3.36AU]), only 1 trial led to the innermost planet
having a semi-major axis less than 1AU, so statistically it is an outlier.
# Solar System Simulation Code:

The code is written in Python 3 and uses the Nbody integrator provided in class (Rauch & Hamilton 2002)(link:https://www.astro.umd.edu/~rauch/HNBody/). The code
first writes an hnb file (named inputs.hnb) to give to the Nbody integrator. The mass of the star is set
to the mass of our sun and the integrator is set to Runge-Katta which is ideal for close encounters between planets.

The number of planets can be adjusted, but are assumed to all be heavy weights. With the given
number of planets the code then draws random numbers for the masses and distances of the planets using
the same distributions given in the article. The initial eccentricities, inclinations, long ascending nodes, and
mean anomalies are zero. The argument of periapsis is randomly selected from the range 0◦ to 360◦.

Once the file is written a function is called to run hnbody.exe with the written file. Once the integration
is done it writes files containing the orbital elements for each planet. Another function can be called to read
in the elements for each planet. With these elements, the orbits of all the bodies can either be plotted
simultaneously or the orbital elements (specifically the semi-major axis, eccentricity, and inclination) of each
body can be plotted against time. There is a function that animates the solar system as well.
The last function called runNsimulations runs the simulations N times (which is given as an argument).
It takes the average of the initial average mass of system, the final average mass of system, the initial minimum 
semi-major axis, the final minimum semi-major axis, and the final average eccentricity of the planets and returns them.

The code is also capable of returning the number of planets ejected and thus the number of planets
remaining. It does this by looking at when the eccentricity of a planet becomes larger than 1 (since that
would mean the orbit is now unbound). Sometimes the eccentricity of a planet goes over 1 then below it
again in a short period of time, so it hasn’t left the system. Therefore, in order to accurately know when a
planet has been ejected, its eccentricity is averaged over a time period. For example, if the averaging time
is 100 yrs and the simulation ran for 1 Myr then the eccentricity is averaged from 999,900 yrs to 1 Myr. If
the eccentricity during that time period is over 1 the planet is considered to have left the system.

# References:

[1] Adams, F. C., & Laughlin, G. (2003). Migration and dynamical relaxation in crowded systems of giant
planets. Icarus, 163(2), 290–306. doi:10.1016/S0019-1035(03)00081-2

[2] Rasio, F. A., & Ford, E. B. (1996). Dynamical Instabilities and the Formation of Extrasolar Planetary
Systems. Science, 274(5289), 954–956. doi:10.1126/science.274.5289.954

[3] Juric, M., & Tremaine, S. (2008, January). Dynamical Relaxation by Planet-Planet Interactions as the
Origin of Exoplanet Eccentricity Distribution. In D. Fischer, F. A. Rasio, S. E. Thorsett, & A. Wolszczan
(Eds.), Extreme Solar Systems (p. 295).

[4] Pater, I. D., & Lissauer, J. J. (2016). Planetary sciences (2nd ed.). Cambridge University Press.

[5] Weidenschilling, S. J., & Marzari, F. (1996). Gravitational scattering as a possible origin for giant planets
at small stellar distances. Nature, 384(6610), 619–621. doi:10.1038/384619a0
