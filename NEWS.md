## genridge 0.7.1 (2024-11-07)

o Added links to gentalk.pdf
o Improved README and fixed some examples

## genridge 0.7.0 (2023-07-31)

o Converted the package to roxygen2 documentation, correcting some infelicities with S3 methods
o Added an extended README example.
o fix link to genridge paper PDF

## genridge 0.6.8 (2023-07-11)
o a maintenance release, correcting a problem with an \alias{} in one .Rd file

## genridge 0.6.7 (2020-01-07)
o Remove references to ElemStatLearn, now defunct

## genridge 0.6-6 (2017-10-01)
o bump package version, preparing to move to Github

## genridge 0.6-5 (2014-11-24)
o Use rgl:: in calls to rgl package for CRAN
o Fixed long lines in .Rd files

## genridge 0.6-4 (2012-2-10)
o Added contourf()

## genridge 0.6-3 (2011-12-28)
o Minor tweaks to plot3d.ridge
o Added biplot.ridge method for symmetry, showing the PCA vectors in variable space
o Added examples to ridge.Rd using data(manpower, package="bestglm")
o Added Manpower data

## genridge 0.6-2 (2011-12-19)
o Made ridge() an S3 generic with formula and default methods
o Fixed labels problems in biplot.pcaridge and plot3d.ridge
o Added Detroit data with examples showing most methods

## genridge 0.6-1 (2011-12-06)
o Fixed order of classes in pca.ridge
o Made plot3d an S3 generic; added plot3d.pcaridge and plot.pcaridge
o Added biplot.pcaridge() to draw variable vectors on a plot.pcaridge() plot
o Added options to precision.* for greater flexibility
o Added GCV computation to ridge()
o Added which.lambda to select ellipses in plot.ridge(), plot3d.ridge() etc.
o Added labels argument to plot.ridge() and fixed pairs.ridge for this.

## genridge 0.6-0 (2011-12-04)
o Added plot3d.ridge() using rgl, completing the main stages in visualization methods for ridge-like problems

## genridge 0.5-3 (2011-11-30)
o Added precision.ridge() method for measures of precision and shrinkage
o Added vif.ridge() method, extending the car S3 generic for ridge objects
o Now Depends: car

## genridge 0.5-2 (2011-11-25)
o Added Acetylene data
o Added vcov.ridge() method

## genridge 0.5-1 (2011-08-22)
o Fixed Authors@R in DESCRIPTION

## genridge 0.5-0 (2011-08-22)

o Added NAMESPACE
o Initial package version, released to R-Forge
o Added pca.ridge
