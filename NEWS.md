# cvmaPLFAM 0.1.1

-   Added a `NEWS.md` file to track changes to the package.

## Major changes

-   Replaced the previous `predRisk` function with the new `cvpredRisk` function, focusing solely on the implementation of the CVMA method. The outputs of other methods have been removed.
-   Extended the `modelspec` function to support nested model structures.
-   Adjusted the order of functional principal component analysis and cross-validation steps in the CVMA process.
-   Made the previously internal function `plam.fit` visible to users.

## Minor changes

-   Updated documentation.
