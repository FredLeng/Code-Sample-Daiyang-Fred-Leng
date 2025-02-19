# Code-Sample-Daiyang-Fred-Leng

This code is a homework assignment for my Nonparametric Smoothing class. It is a regression fit using Splines to the "lidar" data from the package "SemiPar". As you can see, the plot of the data (first plot) shows local behavior that can and ideally should be modeled with a Spline. A truncated linear polynomial basis is used, although a B-Spline basis is ideal, and there are several packages for implementing it. However, there are no compmutational singularity here and the B-Spline was not introduced at the time of this assignment. The Residual Sum of Squares and Cross Validation Residual Sum of Squares for different values of the tuning parameter are calculated and plotted (second plot). The value of tuning parameter with minimum Cross Validation Residual Sum of Squares is then used in the final regression fit of the data (third plot). 
