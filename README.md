# OMP-Structure-Sequence

Here are the source codes for the article [Inward-facing glycine residues create sharp turns in β-barrel membrane proteins](https://www.sciencedirect.com/science/article/pii/S0005273621001127).

## Details

To investigate how the shape of OmpX β-barrels changes with an increasing number of β-strands, we tracked the shape of OmpX variants during equilibrium simulations. The calculated shapes are consistent with our initial observations; the variants with 8, 10, 12 and 14 β-strands maintain a high shape value while that of the 16-stranded variant drops significantly. The shapes of OmpX variants with 10, 12 and 14 strands are about 0.8 as calculated from 3-μs equilibrium simulations, suggesting that they are close to a circle. The shapes of WT OmpX (8-stranded) are relatively low with a value of about 0.6. However, the shapes of variants with 16 strands are about 0.3, demonstrating that these OmpX variants adopt a nearly flat shape.

<p align="center">
  <img src="https://github.com/zzhang624/OMP-Structure-Sequence/blob/main/figs/f1_new.png">
</p>
Fig. Equilibrium simulations of OmpX variants with different strand numbers. A. Snapshots show five OmpX variants built with different strand numbers in cartoon (side view) and surface (top view) representations colored red (N-terminus) to blue (C-terminus). B. Snapshots show the structure of OmpX variants after 3-μs equilibrium simulations. C. Representative snapshots (bottom view) of two OmpX variants at the end of 3-μs equilibrium simulations. Semi-major (a) and semi-minor (b) axes of OmpX β-barrels are shown. D. The shape of OmpX variants over time in the two equilibrium simulations (blue and orange, respectively). Initial shapes are marked by red dashes. E. Shapes of OmpX variants vs. strand numbers. The data points are calculated from the last 300-ns of each 3-μs simulation.


<p align="center">
  <img src="https://github.com/zzhang624/OMP-Structure-Sequence/blob/main/figs/f3_new.png">
</p>

Fig. Location of glycines plays a role in determining the overall shape of OmpX β-barrels. A. Side and bottom views of the 16-stranded OmpX variant β-barrels. Two sets of glycines parallel to the membrane normal, labeled I and II, are shown in a space-filling representation in the side view and marked with circles in the bottom view. B. Schematic definition of the 3-Cα angle. Cα atoms are identified as nearby by virtue of the adjacent backbone hydrogen bonds (green dashed lines) between their β-strands. The 3-Cα angle (θ) was then defined as the membrane-plane projection of the angle between three nearby Cα atoms and was assigned to the central amino acid. All 3-Cα angles are calculated from the inward-facing residues. C. Averaged 3-Cα angles of all amino acids calculated from MD simulations of OmpX variants. Error bars denote ±1 standard deviation. Glycine is shown with two data points reflecting two apparent populations. D. The 3-Cα angle distribution of inward-facing glycines (GLY) in simulations.


<p align="center">
  <img src="https://github.com/zzhang624/OMP-Structure-Sequence/blob/main/figs/f5.png">
</p>

Fig. 3-Cα angles for inward-facing residues. A. Average 3-Cα angles of different amino acids facing inside from 119 β-barrel membrane-protein PDBs. Error bars denote ±1 standard deviation. Glycine is shown with two data points reflecting two apparent populations. B. The distribution of 3-Cα angles for inward-facing glycines. C. Top and side view of β-barrels of OprN (PDB ID: 5AZP), Omp33 (PDB ID: 6GIE), and Ail (PDB ID: 3QRA). Glycines at the sharp turns are highlighted.
