## Test environments
* local Windows 7 install, R 3.4.1
* win-builder R Under development (unstable) (2017-09-12 r73242), 3.4.1 (2017-06-30)
* R-Forge R version 3.4.1 Patched (2017-09-11 r73238)

## R CMD check results
There were no ERRORs or WARNINGS.  There was one NOTE:

Examples with CPU or elapsed time > 10s
         user system elapsed
Detroit 11.89   0.12   12.42

First, I tried wrapping more than half the example in \donttest{},
but it didn't have much impact.
Moreover,  I can't reproduce it in my local session.
When I run the **all** the code in the example, I get:

   user  system elapsed 
   1.66    0.25    2.00 


## Comments
This is a minor release, ...
