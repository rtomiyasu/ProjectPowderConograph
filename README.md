[to Japanese](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/README-jp.md)

# Instructions for Conograph CUI program
[Version 0.99](https://github.com/rtomiyasu/ProjectPowderConograph/tree/main/Conograph1_0_00_win) is the initial release of the core program of the Conograph software which executes powder auto-indexing by the method introduced in [[2](#References)], [[3](#References)] for the first time.

<img alt="outline" src="https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/7ee286db-aee6-40a3-8613-66fa6b806a69" width="50%">
![outline](https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/7ee286db-aee6-40a3-8613-66fa6b806a69)
```
Figure 1: three main stages of powder auto-indexing
```

Conograph executes extensive powder auto-indexing searches in comparatively short time frames for powder diffraction patterns including ones from spallation neutrons. After the peak search is finished, the powder auto-indexing is completed by carrying out the three main stages that are shown in Figure 1. An error-stable Bravais lattice determination is necessary, because Conograph uses an enumeration algorithm that commonly works for any Bravais lattice, space group and type of systematic absences, and observational errors in peak positions are propagated to the unit-cell parameters that are obtained.

When the regular search option is selected, Conograph completes all of the procedures in about 10 minutes using the default parameters prepared for Conograph. In the quick search option, the time will not be more than 5 minutes. (The times were obtained using an Intel® Core™ i7 Processor (3.2 GHz, 8 threads). If a less powerful processor is used, then the search options will take more time.)

If you have not used Conograph before, then it is best to choose the regular search option, because this allows you to use the default searching parameters without changes even in difficult cases. The quick search option without modifications to the parameters also works well if the unit-cell has a small volume or high symmetry. The differences between the two searches is explained [here](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/doc/TipsForPowderAutoIndexing.md#What_are_the_differences_between_the_quick_and_regular_search_options).

## NEWS
### 2016/9/7
- The output format for base-centered monoclinic cells was corrected.

## FAQ
- [How to use the Conograph Command User Interface (CUI)](#How_to_use_the_Conograph_Command_User_Interface)
    - [Starting Conograph](#Starting_Conograph)
    - [Refinement of the unit-cell parameters and the zero-point shift after the execution of a powder auto-indexing search](#Refinement_of_the_unit_cell_parameters_and_the_zero_point_shift_after_the_execution_of_a_powder_auto_indexing_search)
- [Tips for peak searching](#Tips_for_peak_searching)
- [Tips for powder auto-indexing using Conograph](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/doc/TipsForPowderAutoIndexing.md) (link to another web page)
    - A table of parameters in the *.inp.xml file
    - What are the differences between the quick and regular search options?
    - Which parameters should be modified when the results are not satisfactory?
    - The powder auto-indexing outputs many unit-cells. How can I find the correct one? (Properties of figures of merit)
    - Uniqueness of solutions in powder auto-indexing
- [Differences between CUI and GUI](#Differences_between_CUI_and_GUI)

### How_to_use_the_Conograph_Command_User_Interface
#### Starting_Conograph
1. The Conograph program requires the following three input files. (Examples can be found in the sample folder.)
    - A *.inp.xml file that includes information about the input parameters ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/Conograph1_0_00_win/sample/sample5/Cimetidine-SR.inp.xml))
    - A cntl.inp.xml that includes the names of the *.inp.xml file and the output file ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/Conograph1_0_00_win/sample/sample5/cntl.inp.xml))
    - An IGOR text file that includes a powder diffraction pattern (X and Y coordinates and the error in Y) and the following information about peaks ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/Conograph1_0_00_win/sample/sample5/Cimetidine-SR_pks.histogramIgor). This file is outputted from the [peak search program](https://github.com/rtomiyasu/PeakSearch/tree/main) distributed separately)
        1. peak-positions (2θ, time-of-flight, d)
        1. peak-heights (used only for display)
        1. full-widths-at-half-maximum (used to estimate errors of the peak-positions)
        1. 0/1 flags for the respective peaks for powder auto-indexing and computation of figures of merit
1. Copy one of the folders from the sample folder. Modify the contents of the two xml-files, the name of *.inp.xml file, and the 0/1 flags in the IGOR text file if necessary. If you change the name of the *.inp.xml file, then it will be necessary to modify the contents of the cntl.inp.xml file accordingly.
1. Open a command prompt or terminal window in your operating system. Change the current folder to the folder that contains the modified cntl.inp.xml file.
1. Enter the absolute path to the Conograph.exe file on the command line and execute Conograph.
1. The CUI outputs an XML file containing a list of candidate unit-cells ([Example of a hexagonal case](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2.index.xml)) and awaits input from the user. At the top of the XML output file, unit-cells with the best value for the figures of merit are presented for each Bravais type. Just below the unit-cells, the unit-cell with the best value for the de Wolff figure of merit [[4](#References)] from all of the obtained unit cells is presented.
#### Refinement_of_the_unit_cell_parameters_and_the_zero_point_shift_after_the_execution_of_a_powder_auto_indexing_search
1. A number like 0403003 is associated with every unit-cell in the XML output file. If the user enters one of the numbers on the command line, the the unit-cell parameters and zero-point shift parameter are refined and the program outputs an IGOR text file ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2_lattice(Hexagonal_23.1%2C23.1%2C10.7%2C90%2C90%2C120_70.3).histogramIgor)) that includes the following information:
    1. copy of the contents in the input IGOR text file
    1. positions of computed lines of the input unit-cell parameters
    1. the parameters and the Bravais lattice of the unit-cell
1. Until you terminate Conograph by entering quit, you can enter the number associated with any of the unit cells in the XML file as many times as you want. Before termination, Conograph outputs the following:
    - an XML file including information about refined unit-cell parameters and zero-point shifts ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2.index2.xml))
    - an IGOR text file unifying the IGOR text files output during refinement ([Example](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2_lattices.histogramIgor))

### Tips_for_peak_searching
If you have not used Conograph before, then our advice in peak search is to pick up all the diffraction peaks as uniformly as possible on the basis of their peak-heights. (Such results will be obtained easily by using the [peak search program](https://github.com/rtomiyasu/PeakSearch/tree/main) equipped with Conograph.) The enumeration algorithm of Conograph then determines in a comparatively short time which combinations of peaks provide better solutions. Unless you have any prior information about diffraction peaks, it is better to avoid reducing the quality of input information by selecting peaks or removing overlapped peaks artificiality.

### Differences_between_CUI_and_GUI
As long as the IGOR Pro software is available, powder auto-indexing in the CUI is not difficult. The following are the chief advantages of the GUI.

- A peak search program is incorporated in the GUI.
- At the start-up of the GUI program, recommended values are automatically entered in the text boxes for the input parameters.
- Visual comparison between the powder diffraction pattern and the computed lines of the unit-cells is easier in the GUI.
- In the CUI, solutions are sorted by the de Wolff figure of merit. In the GUI, the other sorting criteria can be chosen anytime after the enumeration stage.
- The CUI executes the various functions as successive procedures. In the GUI, the computational time required for repeated execution may be shorter, because the various functions are available separately and independently. These functions include:
estimation of the zero-point shift Δ2θ by the reflection pair method [[1](#References)] (The CUI program outputs several candidate values on the command line in the initial stage of powder auto-indexing.)
estimation of n required for computation of the de Wolff figure of merit Mn (Some unit cells require n > 20 owing to existence of a dominant zone.)
- In order to use one of the above estimated values, it is necessary to terminate the CUI program once. Regarding a., this limitation of the CUI program is partly because use of Δ2θ = 0 is sufficient for many cases. (For example, Conograph has obtained successful results even in the case where the correct Δ2θ is 0.195°.)

Figure 2 is a flowchart of the CUI:

![flowchart](https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/054785a7-46e1-41d3-a7d7-fac3bce1dfa2)

## References
1. C. Dong, F. Wu, H. Chen,<br>Correction of zero shift in powder diffraction patterns using the reflection-pair method, J. Appl. Cryst., 32, pp. 850-853 (1999).
1. R. Oishi-Tomiyasu,<br>Distribution rules of crystallographic systematic absences on the Conway topograph and their application to powder auto-indexing, preprint.
1. R. Oishi-Tomiyasu,<br>Rapid Bravais-lattice determination algorithm for lattice constants containing large observation errors, Acta Cryst. A, 68, pp. 525-535 (2012).
1. P. M. de Wolff,<br>A simplified criterion for the reliability of a powder pattern indexing, J. Appl. Cryst., 1, pp. 108-113 (1968).
