# Random Matrix Theories

This repository contains computational explorations of **Random Matrix Theory (RMT)** and **quantum chaos**, through both theoretical models and empirical nuclear spectra.

The project studies how complex quantum systems exhibit spectral statistics predicted by random matrix ensembles.
---

## Projects

### 1. Coupled Kicked Top

This project implements simulations of the **bipartite coupled kicked top**, a standard model used to study quantum chaos.

The model consists of periodically kicked angular momentum systems whose dynamics transition from regular to chaotic depending on system parameters. The implementation follows work by **Arul Lakshminarayan and collaborators** on coupled kicked tops and the study of entanglement saturation in classically chaotic systems.

The simulations analyze the **Floquet operator spectrum** and compute statistical observables such as:

- Nearest neighbour spacing distribution
- Spectral density

These statistics are compared with Wishart-Lagaurre random matrix predictions (like Marchenko-Pastur for density of eigenvalues) to quantify saturation value of entanglement as the coupled kicked top becomes chaotic as a result of changing parameters like kick strength (or torsion)

---

### 2. Nuclear Chaos Analysis

This project analyzes **empirical nuclear energy level spectra** using Random Matrix Theory.

Energy levels belonging to the same **spin-parity sequence \(J^\pi\)** are extracted from experimental databases (specifically NATIONAL NUCLEAR DATA CENTRE (NNDC)) and analyzed using standard RMT techniques.

The workflow includes:

1. Cleaning and sorting nuclear energy levels  
2. Spectral unfolding  
3. Computing normalised level spacings  
4. Comparing statistics with GOE and Poisson predictions  

Computed observables include:

- Nearest Neighbour Spacing Distribution (NNSD)
- Spacing Cumulative Distribution Function (CDF)
- Spectral Rigidity \( \Delta_3(L) \)
- Spectral Form Factor (SFF)

These analyses test whether nuclear spectra exhibit **chaotic level statistics** consistent with Random Matrix Theory.

---

## References

•⁠  ⁠M. L. Mehta — Random Matrices  
•⁠  ⁠F. Haake — Quantum Signatures of Chaos  
•⁠  ⁠Bohigas, Giannoni & Schmit (1984)  
•⁠  ⁠Lakshminarayan et al. — studies on kicked tops and quantum chaos  

---

## Author

*Gurbani Mahajan*
