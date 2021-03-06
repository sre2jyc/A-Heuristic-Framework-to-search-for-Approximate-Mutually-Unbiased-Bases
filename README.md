### A Heuristic Framework to search for Approximate Mutually Unbiased Bases, CSCML 2022

This paper is accepted as a regular paper in **The 6th International Symposium on Cyber Security, Cryptology and Machine Learning (CSCML 2022)**.

**Paper Link :**

## Authors
* Sreejit Chaudhury, Jadavpur University.

* Ajeet Kumar, Applied Statistical Unit, Indian Statistical Institute.

* Subhamoy Maitra, Applied Statistical Unit, Indian Statistical Institute.

* Somjit Roy, Department of Statistics, University of Calcutta.

* Sourav Sen Gupta, Nanyang Technological University.


_______________________________________________________________________________________________________________________________________

## Abstract
Mutually Unbiased Bases (MUBs) have many applications including Quantum Key Distribution (QKD), which is one of the most well studied domains in quantum cryptology. Construction of MUBs is a challenging problem that received serious attention for several decades, from theoretical interests too and number of questions are unsolved till date. As set or required number of MUBs may not always be available for different composite dimensions, Approximate MUBs (AMUBs) received serious attention in literature. In this paper, we explain a heuristic to obtain AMUBs with significantly good parameters. Given a non-prime dimension d, we note the closest prime d' > d and form d'+1 MUBs through well known techniques. Then our proposed idea is (i) to apply basis reduction techniques that are well studied in Machine Learning literature in obtaining initial solutions, and finally (ii) to exploit steepest ascent kind of search to obtain further improved results. The efficacy of our technique is shown through construction of AMUBs in dimensions d = 6, 10, 46 from d' = 7, 11 and 47 respectively. Our technique provides a novel framework in construction of AMUBs that can be refined in a case-specific manner. From a more generic view, this approach considers approximately solving a challenging (where efficient deterministic algorithms are not known) mathematical problem in discrete domain through state-of-the-art heuristic ideas.


_______________________________________________________________________________________________________________________________________

## Approach
For Merged Dimension Reduction Approach -
![Merged Approach](/images/MergedSchematic.png)

For Non-Merged Dimension Reduction Approach -
![Non-Merged Approach](/images/NonMergedSchematic.png)


_______________________________________________________________________________________________________________________________________

## Results
Table 1 : Numerical Results For Merging Technique : Algorithm 1.
![Table1](/images/Table1.png)

Table 2 : Numerical Results For Heuristic Search (Algorithm 4) on the bases obtained from Merged Technique wrt ASD (D^2).
![Table2](/images/Table2.png)

Table 3 : Numerical Results For Heuristic Search (Algorithm 4) on the bases obtained from Merged Technique wrt Drift Measure (S).
![Table3](/images/Table3.png)

Table 4 : Numerical Results For Non-Merging Technique : Algorithm 2.
![Table4](/images/Table4.png)

Table 5 : Numerical Results For Heuristic Search (Algorithm 4) on the bases obtained from Non-Merged Technique wrt ASD (D^2).
![Table5](/images/Table5.png)

Table 6 : Numerical Results For Heuristic Search (Algorithm 4) on the bases obtained from Non-Merged Technique wrt Drift Measure (S).
![Table6](/images/Table6.png)

_______________________________________________________________________________________________________________________________________
## Code Guide

```
+-- Heuristics
|   +-- Merged_4
|   +-- Merged_5
|   +-- NonMerged_4
|   +-- NonMerged_5
|   +-- NumberOfIter
|   +-- RandomBases
|   +-- SecondMeasure
|   +-- .DS_Store
+-- Merged_4sets
|   +-- dim6
    |   +-- dim6.py
|   +-- dim10
    |   +-- dim10.py
|   +-- dim46
    |   +-- dim46.py
+-- Merged_5sets
|   +-- dim6
    |   +-- dim6.py
|   +-- dim10
    |   +-- dim10.py
|   +-- dim46
    |   +-- dim46.py
+-- NonMerged_4sets
|   +-- dim6
    |   +-- dim6.py
|   +-- dim10
    |   +-- dim10.py
|   +-- dim46
    |   +-- dim46.py
+-- NonMerged_5sets
|   +-- dim6
    |   +-- dim6.py
|   +-- dim10
    |   +-- dim10.py
|   +-- dim46
    |   +-- dim46.py
+-- images 
|   +-- MergedSchematic.png
|   +-- NonMergedSchematic.png
+-- README.md
````


_______________________________________________________________________________________________________________________________________

