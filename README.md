# scConform: uncertainty quantification for cell type annotation using conformal inference 

## Install from github
The package can then be installed via:

```
if (BiocManager::version() >= "3.23") {
    BiocManager::install("scConform")
} else {
    spdl::info(
        "'scConform' requires Bioconductor version 3.23 or later, ",
        "installing development version from Github"
    )
    devtools::install_github("ccb-hms/scConform")
}
```

Load the library.
```
library(scConform)
```
