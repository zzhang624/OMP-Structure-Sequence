# OMP-Structure-Sequence

Here are the source codes for the article [Inward-facing glycine residues create sharp turns in β-barrel membrane proteins](https://www.sciencedirect.com/science/article/pii/S0005273621001127).

## Shape calculation (projection.tcl & get_b_over_a.py)
For each OmpX variant, an ellipse is fit to the x and y coordinates of the β-barrel Cα atoms, which are flattened onto the z=0 plane, using the cv2 python library and cv2.fitEllipse function. We extract the semi-major (a) and semi-minor (b) axes from these ellipses and define the ratio of b to a as the shape of the fitted β-barrel. Thus, the shape value will fall between 0 and 1, with a perfect circle giving a shape value of 1 and a flattened ellipse giving a shape value close to zero.

## 3-Ca angle 
We downloaded all β-barrel membrane-protein structures analyzed in this paper from the Orientation of Proteins in Membranes (OPM) server. To sort and filter these structures, we first searched for individual structure information using OPM's API (Application Programming Interface) and selected only those structures that included a β-strand number. Then, we further refined the list by removing any apparent duplicates. Specifically, we defined two structures as duplicates if they were the same protein from the same species but in different conformations. If the same protein had structures from two different species, they were considered two separate structures. In total, we analyzed 119 different OMP structures. Next, we identified those residues in the β-barrel of the protein by selecting protein residues between dummy atoms that mark the transmembrane region of the membrane proteins from the OPM server. Then, we manually inspected each structure to exclude any protein residues that are inside the barrel but are not part of it. The barrel-only regions of OMPs were further processed using VMD's psfgen plugin to add hydrogens. To determine the orientation of each residue's side chain, we first imposed two different vectors on the surface area of the barrel: one from the center of the barrel to the backbone of the residue and the other from the center of the barrel to the side chain of the residue. If the backbone vector was shorter than the side-chain vector, then the residue was considered outward-facing; conversely, if the side-chain vector was shorter than the backbone vector, the residue was considered inward-facing. For glycine residues, to define the backbone and the side chain accurately, the backbone selection was changed to the HA2 atom and the side-chain selection to the HA1 atom (CHARMM force-field names). We also utilized the expectation that residues in β-strands alternate facing outward and inward to confirm the in/outward recognition results. 

## Results

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
