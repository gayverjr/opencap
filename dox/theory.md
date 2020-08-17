Theory {#theory}
================

Resonances and Non-Hermitian Quantum Mechanics
----------------------------------------------
Electronic resonances are metastable electronic states with finite lifetimes embedded in the
ionization/detachment continuum. Common examples include temporary anions formed by
electron attachment, and core-excited and core-ionized states which can undergo Auger decay or similar 
relaxation pathways. These states are not part of the usual \f$L^2\f$ Hilbert space of
square integrable functions, and instead belong to the continuous spectrum of the electronic Hamiltonian. 
Theoretical description of resonances is generally not possible by means of conventional 
bound-state quantum chemistry methods, and special techniques are required to obtain accurate 
energies and lifetimes.

Non-Hermitian quantum mechanics (NHQM) techniques provide an attractive approach
that enables adaptation of existing quantum chemistry methodologies to treat metastable electronic
states. In NHQM formalisms, a resonance appears as a single square-integrable
eigenstate of a non-Hermitian Hamiltonian, associated with a with a complex eigenvalue: 

\f$ E=E_{res}-i\Gamma/2 \f$
 
The real part of the energy \f$ (E_{res}) \f$ is the resonance position. The imaginary part \f$ (\Gamma/2) \f$ 
is the half-width, which is inversely proportional to the lifetime of the state[1].

Complex Absorbing Potential
---------------------------
Complex absorbing potentials (CAPs) are imaginary potentials added to the Hamiltonian, and
they are routinely used for evaluation of resonance parameters. In this context, %CAPs 
transform a resonance into a single square integrable state, rendering it accessible by 
means of standard bound-state techniques. To this end, the electronic Hamiltonian is 
augmented with an imaginary potential:

\f$ H^{CAP}=H-i\eta W \f$

where \f$ \eta \f$ is the %CAP strength parameter, and W is a real potential which vanishes in the
vicinity of the molecular system and grows with distance[2].

Since the %CAP-augmented Hamiltonian depends on the strength of the %CAP, a choice
has to be made on the optimal value of \f$ \eta \f$ that provides best estimates of the resonance
position and width. In a complete one-electron basis, the exact resonance position and
width are obtained in the limit of an infinitesimally weak %CAP \f$ (\eta \rightarrow 0^+) \f$. In practice
when finite bases are used, an optimal %CAP strength \f$ \eta_{opt} \f$ is found by locating a stationary
point on the eigenvalue trajectory E(\f$ \eta \f$). A commonly used criterion is 
the minimum of the logarithmic velocity (\f$|\eta\frac{dE}{d\eta}|\rightarrow min \f$)[2].


Perturbative or "Projected" CAP
----------------------------------------
There are multiple strategies for how to incorporate %CAPs into an electronic structure calculation. 
The most straightforward implementation is to engage the one-electron %CAP term starting at the 
lowest level of theory (e.g. Hartree-Fock). While conceptually simple, this requires 
modification of electronic structure routines to handle the complex objects. 
Additionally, this approach requires a unique calculation for each \f$ \eta \f$ along the 
eigenvalue trajectory, which can become prohibitively expensive for larger systems or 
dynamical simulation. 

An efficient alternative is to treat the %CAP as a first order perturbation, considering only a 
small subset of the eigenstates of the real Hamiltonian[3]. In this case, the %CAP will be 
introduced in the basis of the reduced subset of states:

\f$ W_{uv}=\langle u | W | v \rangle \f$

where \f$u\f$ and \f$ v \f$ are eigenstates of the real Hamiltonian. Since the %CAP is a 
one-particle operator, these expressions can easily be evaluated using 
the %CAP matrix in atomic orbital basis evaluated separately, the one-electron reduced 
density matrices (\f$ \rho \f$) for each state, and the set of transition density matrices
(\f$ \gamma \f$) between each pair of states that are obtained from the bound-state 
calculation.

\f$ W_{uv}= Tr\left[W^{AO}\gamma^{uv} \right ], u \neq v \f$

\f$ W_{uv}= Tr\left[W^{AO}\rho^{u} \right ], u = v \f$
 
Once %CAP matrix is evaluated the %CAP-augmented Hamiltonian is constructed as follows:

\f$ H^{CAP}=H_0-i\eta W \f$

where \f$ H_0 \f$ is an appropriate zeroth order Hamiltonian obtained from 
the electronic structure calculation, and \f$ W \f$ is the %CAP represented in the subspace. 
Diagonalization of this %CAP-augmented Hamiltonian
yields \f$ \eta \f$-depdendent eigenvalues that are used to extract 
resonance position and width. Importantly, as only a small number of states in considered 
(typically less than 30), finding the eigenvalues of the %CAP-augmented Hamiltonian has negligible cost in 
comparison to the bound-state electronic structure calculation required to get the initial
set of states (u,v,..). Thus although this perturbative or projected approach introduces 
another parameter (number of eigenstates), the overall cost is essentially reduced to that
of a single electronic structure calculation.
 
With the zeroth order Hamiltonian and the %CAP matrix, eigenvalue trajectories can be 
generated by means of simple external scripts, and estimates of resonances positions and 
widths can be obtained from the analysis of the trajectories.

References
-----------
[1] Reinhardt, W. P. Complex Coordinates in the Theory of Atomic and Molecular Structure and Dynamics. *Annu. Rev. Phys. Chem.* **1982**, 33 (1), 223–255.

[2] Riss, U. V.; Meyer, H. D. Calculation of Resonance Energies and Widths Using the Complex Absorbing Potential Method. *J. Phys. B At. Mol. Opt. Phys.* **1993**

[3] Sommerfeld, T.; Santra, R. Efficient Method to Perform CAP/CI Calculations for Temporary Anions. *Int. J. Quantum Chem.* **2001**, 82 (5), 218–226.
