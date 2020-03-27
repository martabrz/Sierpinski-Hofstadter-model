# topological-fractals
Codes for the paper 'Topology in the Sierpiński-Hofstadter problem' (Phys. Rev. B 98, 205116, (2018)). Written in C++ with Armadillo (http://arma.sourceforge.net/) library. Tested on clang version 11.0.0.


### What's here?
There are three main folders: SierpinskiCarpet, SierpinskiGasket and PlottingScripts (complementary scripts in Python3 for publication-quality plots), together with an example of makefile.

```c++
CreateSierpinskiCarpet
```
Creates the Sierpiński carpet as a lattice at n-th iteration based on the Hadamard matrix method.

```c++
CreateSierpinskiTriangle
```
Creates the Sierpiński carpet as a lattice at n-th iteration based on the Pascal triangle mod integer method.

```c++
CreateHamiltonian
```
Using previously created lattice, it builds a nearest-neighbour tight-binding Hamiltonian for a given hopping strength t and flux value (+ separate function OnSiteDisorder to study the effect of disorder).


### Methods
```c++
LocalChern
```
Computes the Chern number in a real-space (see 'Anyons in an exactly solved model and beyond' by A. Kitaev). Checks the quantization as a function of a patch size.

```c++
BottIndex
```
Computes the Bott index (see T. A. Loring and M. B. Hastings, EPL92, 6 (2011)).

```c++
IPR
```
Inverse participation ratio to define whether states are localized or delocalized.

```c++
LocalizationSpatial
```
|\psi|^2 of the wavefunction to study the edge modes.

```c++
LevelStatistics
```
Study topological phase transition as a function of disorder by computing the variances from level statistics.

```c++
EdgeLocality
```
Investigate the eigenstates localization properties with the edge-locality marker.
