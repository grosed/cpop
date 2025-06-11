#  **News**


## cpop 1.0.8

- Added code to ensure variance >= epsilon = 1e-6. Many thanks to Piotr Fryzlewicz for raising this issue.

## cpop 1.0.7

- Added CITATION file
- Removed vinettes. Please refer to doi:10.18637/jss.v109.i07 instead.
- Removed Makevars and Makevar.win (assume the user has c++ >= 14 available)
- Added reference to doi:10.18637/jss.v109.i07 in documentation


## cpop 1.0.6

- Updated vignettes to be consistent with article submitted to the Journal of Statistical Software.

## cpop 1.0.5

- Use message to return messages instead of cat.

## cpop 1.0.4

- Added NEWS.md
- Added guards for setting show and fitted methods as generic to preventing masking
- Removed redundant (non exported) null2na functions
- Changed name of vignettes from article to cpop
- Corrected details of reference in DESCRIPTION file
- Changed titles in documentation to sentence case
- Improved documentation of cpop method
- Changed default value of sd in cpop to be data dependent
- Fixed error in estimate function

##  cpop 1.0.3 

### Updated on CRAN (2022-08-23)

- Added Vignettes
- Renamed simulate function to simchangeslope
- Fixed minor error in wave_number data set documentation
- Added fix for issue regarding nearly singular design matrices
- Added extra documentation regarding default values for beta in cpop
- Renamed sigma parameter to sd in simulate 
- Updated documentation and examples for cpop.crops method
- Fixed some broken URLs


##  cpop 1.0.0 

### Initial release to CRAN  (2022-01-11)
