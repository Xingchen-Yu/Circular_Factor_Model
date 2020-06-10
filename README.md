# Circular_Factor_Model

## Code
Circular_Factor_Model.R takes 4 arguments which are iterations, burn-in period of mcmc, number of cores used, the number of the U.S
House of representatives data. It is recommended to run on a server with many cores.

For exammple, to run 100000 iterations of which 80000 are burn-in using 12 cores for the 112th U.S House of representative data,
one should excute the following commands.

Rscript ./main_script/Circular_Factor_Model.R 100000 80000 12 112 &




