# MavericK

Here are the raw files that make up the program [MavericK](www.bobverity.com/maverick). These files are distributed under [MIT license](https://opensource.org/licenses/MIT), which means you can basically do whatever you want with them as long as you include the original copyright and license notice in any copy of the software/source code.

The Makefile contains all the commands needed to compile the program on Unix-like systems (i.e. Mac and Linux), meaning you should only need to implement the command "make" on the command line to compile the program (this has been found to work using GCC version 6.1 and above). On a Windows machine it is recommended to load all .cpp and .h files into Visual Studio where they can be compiled easily.

Separate branches exist for each released version of the software. Select the branch corresponding to the version you are interested in before downloading files, or alternatively you can see all released versions in the "releases" tab.


## Version History

12.12.2017  **Version1.0.5 released**<br>
11.12.2017  Allow alternative data file format with multiple columns per locus.

07.12.2016  **Version1.0.4 released**<br>
30.11.2016  Increased speed under admixture model via look-up table.

29.06.2016  **Version1.0.3 released**<br>
29.06.2016  Improved cross-platform compatibility. Fixed clock time, which previously gave incorrect runtime estimates.

22.06.2016  **Version1.0.2 released**<br>
21.06.2016  Added Metropolis-Hastings step to improve mixing under admixture model

20.05.2016  **Version1.0.1 released**<br>
20.05.2016  Fixed bug causing Q-matrix probabilities to sum to >1

20.05.2016  **Version1.0.0 released**
