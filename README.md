# GalaxySim
Simulation of the tidal interaction of two galaxies for different initial conditions such as
* Light Intruder
* Heavy Intruder
* Retrograde Passage
* Direct Passage

based off of the 1972 Toomre and Toomre paper
https://pubs.giss.nasa.gov/abs/to03000u.html


# Galaxy Mergers Simulation

## Introduction
Galaxy mergers play a crucial role in the evolution of galaxies, leading to phenomena such as tidal tails and galactic bridges. This project aims to simulate galaxy mergers using numerical methods, following a path similar to the work of the Toomre Brothers and subsequent research.

## Background
In recent years, scientists have utilized numerical methods to accurately model galaxy interactions. The Toomre Brothers' simulation in FORTRAN demonstrated the formation of tidal tails and galactic bridges, laying the foundation for later studies. Subsequent research, such as Bournaud, Duc, & Emsellem's work in 2008, utilized advanced supercomputers to model millions of particles per galaxy, revealing further insights into the merger process.

## Methodology
Our simulation involves two galaxies composed of evenly distributed stars at multiple radii. Using Matlab's ode45 function, we solve Newton's Law of Gravitation to model the gravitational interaction between the galaxies. The simulation explores various scenarios, including retrograde and direct passages of intruder galaxies, with different mass distributions.

## Results
### Retrograde Passage
Equal mass galaxies with opposite spin directions resulted in the intruder galaxy accreting stars from the parent, forming a tidal tail. The parent galaxy also experienced alterations in radius.

### Direct Passage
Equal mass galaxies with the same spin direction demonstrated symmetric disruptions, consistent with Newton's law of Gravitation.

### Lighter Intruder
A lighter intruder in direct passage disrupted both galaxies, with more pronounced effects on the lighter galaxy.

### Heavier Intruder
A heavier intruder in direct passage disrupted both galaxies, with significant disruptions observed in both.

## Conclusion
This project recreates aspects of the Toomre Brothers' work, demonstrating the application of numerical methods in understanding galaxy mergers. While not capturing all complexities, the simulation confirms the formation of tidal tails and bridges observed in previous studies, providing valuable insights into galaxy evolution.

## References
1. Toomre Brothers Simulation: [Link](http://www.spaceref.com/news/viewpr.html?pid=8197)
2. Bournaud, Duc, & Emsellem Paper: [Link](https://academic.oup.com/mnrasl/article/389/1/L8/996886)
