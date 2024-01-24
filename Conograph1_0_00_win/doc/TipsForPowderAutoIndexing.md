[to Japanese](https://github.com/arimuratak/Conograph-powder-auto-indexing-/blob/main/Conograph1_0_00_win/doc/TipsForPowderAutoIndexing_JP.md)

# Tips for powder auto-indexing using Conograph
- [A table of parameters in the *.inp.xml file](#Table_1_Parameter)
- [What are the differences between the quick and regular search options?](#What_are_the_differences_between_the_quick_and_regular_search_options)
- [Which parameters should be modified when results are not satisfactory?](#Which_parameters_should_be_modified_when_results_are_not_satisfactory)
    1. For a more careful search:
        - SearchLevel (XML),
        - CriticalValueForLinearSum (XML),
        - MaxNumberOfPeaks (XML).
    1. For improvement of the efficacy of the figures of merit:
        - ZeroPointShiftParameter (XML),
        - MaxNumberOfPeaksForFOM (XML),
        - 0/1 flags on use of the respective peaks for powder auto-indexing and computation of figures of merit (IGOR text).
    1. For speeding up the process:
        - NumberOfThreadsToUse (XML),
        - MinPrimitiveUnitCellVolume, MaxPrimitiveUnitCellVolume, MinDistanceBetweenLatticePoints (XML),
        - OutputCubicF, OutputCubicI, OutputCubicP, OutputHexagonal, OutputRhombohedral, OutputTetragonalI, OutputTetragonalP, - - OutputOrthorhombicF, OutputOrthorhombicI, OutputOrthorhombicC, OutputOrthorhombicP, OutputMonoclinicB, OutputMonoclinicP (XML),
        - SearchLevel & MaxNumberOfLatticeCandidates (XML).
- [The powder auto-indexing outputs many unit-cells. How can I find the correct one? (Properties of figures of merit)](#Properties_of_figures_of_merit)
- [Uniqueness of solutions in powder auto-indexing](#Uniqueness_of_solutions_in_powder_auto_indexing)

Conograph users are counseled to start powder auto-indexing with the recommended values from the following table for all of the parameters, except for the ones with an empty, gray box. When AUTO is entered, the recommended value for the parameter is computed using the powder diffraction pattern and is presented on the command line.

## Table_1_Parameter
| XML tag | Explanation |Recommended value |
|----------|---------------|:-------------------:|
|Parameters for diffractometer| | |
|IsAngleDispersion　　　　　| 0: time-of-flight, 1: angle dispersion   |                   |
| ZeroPointShiftParameter　| (For angle dispersion) Zero-point shift parameter Δ2θ (°)  | 0 |
| WaveLength 　　　　　　　| (For angle dispersion) Wave length (Å).　|                   |
| ConversionParameters    | (For time-of-flight) Conversion parameters,?i.e., coefficients of a polynomial of any degree starting from the constant term   |                   |
|Parameters for enumeration stage| | |
|SearchLevel|0: quick search (for unit cells with small volume or high symmetry),<br>1: regular search (for all cases including difficult ones).|    |
| MaxNumberOfPeaks | Number of peaks used for calculation | AUTO |
| CriticalValueForLinearSum | The criterion c to judge if linear sums of q-values (including Ito's equation) equals zero (i.e.｜ΣiQi｜<= cErr[ΣiQi]) | 1.0 |
| MinPrimitiveUnitCellVolume,<br>MaxPrimitiveUnitCellVolume | Lower and upper threshold for the volume of primitive unit-cell (Å3) | AUTO |
| MaxNumberOfTwoDimTopographs | Upper threshold for the number of quadruples (q1,q2,q3,q4) taken from selected topographs | AUTO |
| MMaxNumberOfLatticeCandidates | Upper threshold for the number of (triclinic) solutions | AUTO |
| Parameters for Bravais lattice determination |  |  |
| OutputCubicF, OutputCubicI,<br>OutputCubicP, OutputHexagonal,<br>OutputRhombohedral, OutputTetragonalI,<br>OutputTetragonalP, OutputOrthorhombicF,<br>OutputOrthorhombicI, OutputOrthorhombicC,<br>OutputOrthorhombicP, OutputMonoclinicB,<br>OutputMonoclinicP, OutputTriclinic | Output unit cells belonging to the Bravais type? (0:No, 1:Yes) | 1 |
| Parameters for sorting stage |  |  |
| Resolution | Relative resolution of d* = 1/d-values<br>(used for detection of almost identical unit cells) | 0.03 |
| MaxNumberOfPeaksForFOM | Number n of peaks used to calculate figures of Merit. (Even if the number of observed peaks in your diffraction pattern is smaller than 20, it is not necessary to decrease this number.) | 20 |
| MinFOM | Any indexing solution is output if the value of the de Wolff figure of merit Mn is greater than this value. | 3.0 |
| MinNumberOfMillerIndicesInRange,<br>MaxNumberOfMillerIndicesInRange | Lower and upper threshold for the number of distinct computed lines within the range of the first-to-n'th observed lines. | AUTO |
| MinUnitCellEdgeABC | Lower thresholds for the unit cell parameters a, b, c (Å) | 0 |
| MaxUnitCellEdgeABC | Upper thresholds for the unit cell parameters a, b, c (Å) | 1000 |
| Configuration parameters |  | |
| NumberOfThreadsToUse | Number of threads used (Enter MAX to use all threads.) | |
| AxisForRhombohedralSymmetry | Enter Rhombohedral or Hexagonal (rhombohedral axis or hexagonal axis). | |
| AxisForMonoclinicSymmetry | Enter A or B or C (A: a-axis, B : b-axis, C : c-axis). | B |
| ThresholdOnNormM | Unit cells with MnWu less than this value are removed before Bravais lattice determination | 1.9 |
| ThresholdOnRevM | Unit cells with MnRev less than this value are removed before Bravais lattice determination | 1.0 |
| MinDistanceBetweenLatticePoints | If the distance (Å) between two closest points of the crystal lattice is less than this value, the unit cell is removed during the enumeration process. | 2.0 |

## What_are_the_differences_between_the_quick_and_regular_search_options
The two search options were designed to resolve a technical dilemma; namely, whether to give priority to computational speed or to memory usage?. Because these searches are based on the same algorithm, almost identical results will be generated if the input parameters are chosen properly for quick search. (In that case, quick search will be more than twice as fast as regular search.) The following table summarizes the differences between the two searche options:

Table 2: Diferences between the two search options

| | Quick search | Regular search |
|----------|------------|-------------|
| Computational time<br>(Intel® Core™ i7 Processor (3.2 GHz), 8 threads) | < 5 minutes | about 10 minutes |
| Memory usage | A large amount of memory is sometimes required for difficult cases. In order to avoid memory allocation errors, unit cells with a smaller volume are removed first from the list of enumerated solutions, if very many solutions are generated. | Not so much, because unit cells with a better figure of merit are reserved first |
|Modification of input parameters |	Unnecessary in many cases if the unit cell has a small volume or high symmetry. Otherwise it is frequently necessary to modify one of the following parameters: MaxPrimitiveUnitCellVolume or MaxNumberOfLatticeCandidates | Unnecessary in many cases |

## Which_parameters_should_be_modified_when_results_are_not_satisfactory
1. For_a_more_careful_search
    - **SearchLevel: 0 → 1.**<br>Regular search with the recommended values of input parameters will be sufficient for many cases.
    - **CriticalValueForLinearSum: 1 → 1.5.**<br>For some powder diffraction patterns from characteristic X-rays and reactor sources, 1.5 seems to work more efficiently.
    - **MaxNumberOfPeaks: AUTO → a number greater than 48 (50--80).**<br>When AUTO is used, the number entered in this parameter is not greater than 48. Although computational time is magnified by increasing this parameter, more information on peaks is used for powder auto-indexing. As a result, a wider variety of unit cell parameters are formed and checked.
1. For_improvement_of_the_efficacy_of_the_figures_of_merit
    - **ZeroPointShiftParameter: 0 → estimated value of zero point shift.**<br>Large zero-point shift (e.g., Δ2θ) lowers the values of figures of merit largely. If the result of the refinement after powder auto-indexing is insufficient, powder auto-indexing with an estimated value of zero-point shift may be effective. The estimated value is obtained by one of the following method:
        1. the reflection pair method (before powder auto-indexing) [[2](#References)],
        1. refinement of the zero-point shit using unit cell parameters with large MnRev (e.g., >3) and MnSym (e.g., >10) (after powder auto-indexing).0
    - **MaxNumberOfPeaksForFOM: 20 → a number greater than 20.**<br>For this parameter, n = 20 is frequently used, and works well in many cases. However, greater n is required when a domnant zone exists. If a unit cell requiring large n is entered on the command line for refinement, the n and a warning message are displayed on the screen. Information on this n is also included in the XML output files.
    - **0/1 flags on use of the respective peaks for powder auto-indexing and computation of figures of merit: 1 → 0**<br>The sorting of unit cells is sensitive to false peaks caused by impurities. (The enumeration stage, on the other hand, is proved to be insensitive to them.) Therefore, if the MaxNumberOfPeaksForFOM parameter (the number of peaks used for computation of figures of merit) is set to n, the number of false peaks included in the range of the first n peaks should be as small as possible. In order to improve powder auto-indexing results, it is sometimes effective to restart powder auto-indexing from peak search or use 0/1 flags to exclude suspicious peaks. (However, this exclusion may reduce the success rate of enumeration stage to some degree.)
1. For_speeding_up_the_process
    - **NumberOfThreadsToUse** <br>
    A very straightforward method is to increase this parameter.
    - **MinPrimitiveUnitCellVolume, MaxPrimitiveUnitCellVolume, MinDistanceBetweenLatticePoints**<br>
    If you have prior information about these parameters, they will be useful for speeding up of the enumeration stage. (Note that these parameters are concerned with the primitive cell (or primitive lattice) of the crystal, and not with the Bravais lattice.)
    - **SearchLevel: 1 → 0, and MaxNumberOfLatticeCandidates: AUTO → a number greater than 64000 (100000 - 300000)**<br>
    Considering computers with a small amount of memory, the MaxNumberOfLatticeCandidates parameter is set to a number not more than 64000, if AUTO is entered. If the number of unit cells generated is greater than the value in the MaxNumberOfLatticeCandidates parameter, unit cells with a smaller volume are removed first. This is the reason that a modification to either of the paramters MaxPrimitiveUnitCellVolume or MaxNumberOfLatticeCandidates is frequently required for the quick search option. However, for example, MaxNumberOfLatticeCandidates = 300000 did not cause memory allocation errors on a computer with 4 GB memory and the quick search was still at least twice as fast as the regular search. Therefore, by setting the SearchLevel parameter to quick search and increasing the value in the MaxNumberOfLatticeCandidates parameter to a value that is acceptable for your computer, the same result as regular search will be obtained in a shorter time.
    - **OutputCubicF, OutputCubicI, OutputCubicP, OutputHexagonal, OutputRhombohedral, OutputTetragonalI, OutputTetragonalP, OutputOrthorhombicF, OutputOrthorhombicI, OutputOrthorhombicC, OutputOrthorhombicP, OutputMonoclinicB, OutputMonoclinicP: 1 → 0**<br>
    If you have prior information about the Bravais lattice of the crystal, it will also be useful for speeding up of the Bravais lattice determination and the sorting stage.

## Properties_of_figures_of_merit
### The powder auto-indexing outputs many unit-cells. How can I find the correct one?
The most reliable way to distingush the correct unit cell is to compare computed lines of each solution with positions of diffraction peaks visually. By using figures of merit, it is possible to detect more plausible solutions automatically and reduce the number of unit cells to check largely.

The following figures of merit are available at present:

- $`M_n`$ : de Wolff figure of merit [[1](#References)],
- $`{M_n}^{W_u}`$ : improved figure of merit proposed by Wu [[4](#References)],
- $`{M_n}^{Rev}`$ : reversed figure of merit,
- $`{M_n}^{Sym}`$ : symmetric figure of merit,
- $`MN_ε`$ : number of almost identical unit-cells among the ones enumerated,

where n and ε are the values for the MaxNumberOfPeaksForFOM and Resolution parameters, respectively. The values of the first four figures of merit become close to 1, if observed peak positions and computed lines of a fixed unit cell have no correlation. Unit cells with $`M_n > 10`$, $`{M_n}^{W_u}>10`$, $`{M_n}^{Rev}>3`$ or $`{M_n}^{Sym}>30`$ may be considered to be a plausible candidate.
Among the available figures of merit, de Wolff figure of merit $`M_n`$ is superior in efficacy to the others. This is partly because of c, one of the following properties of $`M_n`$:

1. Sensetivity to false peaks caused by impurities
1. Insensitivity to extinct reflections
1. The higher symmetric unit cell always receives a greated value when almost identical solutions exist and belong to different Bravais lattices

Although $`{M_n}^{W_u}`$ also has the properties a and b, it has a marked tendency to select lower symmetic unit cells. $`{M_n}^{Rev}`$ has properties opposite to those of $`M_n`$. In particular, it is insensitive to falase peaks. $`{M_n}^{Sym}`$ is given properties intermediate between $`M_n`$ and $`{M_n}^{Rev}`$.
A unit cell may be considered most likely correct if it received an outstanding value for $`M_n`$ and is also ranked at the top at least among the unit cells of the same Bravais type for the other figures of merit. If every figure of merit ranks a different unit cell at the top in the group of the same Bravais type, various plausible solutions would have been generated. Even so, the correct solution is contained frequently in the list of enumerated unit cells and may be found by checking the list manually.

However, even if the correct solution is included in the list of the enumerated solutions and the observed and computed lines are compared visually, this does not mean that the correct solution will be distinguished, especially in cases where the powder diffraction pattern has poor quality. Regarding this problem, you should also read the paragraphs for [uniquness of solutions in powder auto-indexing](#Uniqueness_of_solutions_in_powder_auto_indexing),

## Uniqueness_of_solutions_in_powder_auto_indexing
Powder auto-indexing solutions are not always uniquely determined from just the information about peak positions. Although existence of more than one solutions occurs in rare instances for low-symmetry unit cells, it is known that cubic, hexagonal, or rhombohedral cells have the same computed lines as lower-symmetry cells [[3](#References)]. This is a consequence of mathematical theorems and is not caused by observational problems. The following figure presents an example of a cubic(F) case:

![ThreeSolutions](https://github.com/arimuratak/Conograph-powder-auto-indexing-/assets/149344913/c8b073f8-c40c-4222-8d86-c1e9244a11f4)
```
Figure 1 : Example of different unit cells with the same computed lines.
```
Even in such a case, all the solutions are considered to be contained in the list of enumerated unit cells if observational problems are not severe because Conograph executes an extensive powder auto-indexing search. The best de Wolff figure of merit will be then given to the unit cell with the highest symmetry. Although this result is consistent with the general belief that higher-symmetry unit cells are more likely to be correct for the powder auto-indexing, it should be noted that this is not always the case. Consideration of systematic absences and the size of unit cells will be helpful for reducing the number of candidates in some cases.

## References
1. P. M. de Wolff,<br>A simplified criterion for the reliability of a powder pattern indexing, J. Appl. Cryst., 1, pp. 108-113 (1968).
1. C. Dong, F. Wu, H. Chen,<br>Correction of zero shift in powder diffraction patterns using the reflection-pair method, J. Appl. Cryst., 32, pp. 850-853 (1999).
1. A. D. Mighell & A. Santoro,<br>Geometrical Ambiguities in the Indexing of Powder Patterns, J. Appl. Cryst., 8, pp. 372-374 (1975).
1. E. Wu,<br>A modification of the de Wolff figure of merit for reliability of powder pattern indexing, J. Appl. Cryst., 21, pp. 530-535 (1988).