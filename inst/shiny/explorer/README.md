# LongOmicsTraj Explorer (Shiny prototype)

Browser-based OTS query builder using the public `lot_glucold_ots` dataset.
No PRISM or other private patient data are included.

## Local launch

Install the package from this branch and optional UI dependencies:

```r
install.packages(c("shiny", "DT", "rsconnect", "remotes"))
remotes::install_github("nzkermani/LongOmicsTraj", ref = "feature/shiny-ots-explorer", dependencies = TRUE)
shiny::runApp("inst/shiny/explorer")
```

## Deployment

Authenticate privately using the instructions in your shinyapps.io account
dashboard. Never commit an rsconnect token or secret.

From the repository root:

```r
rsconnect::deployApp(appDir = "inst/shiny/explorer",
                     appName = "LongOmicsTraj-Explorer")
```

The server must be able to install the GitHub version of LongOmicsTraj;
check the deployment manifest and logs before publishing.

## Scope and caveats

- Queries call the existing `lot_query()` method.
- Assay, estimator, treatment groups and topology labels are read from OTS metadata.
- The templates assign labels to groups in alphabetical order; verify choices.
- Exports include matched OTS rows and the corresponding R query.
- OTS patterns do not establish significance, drug effects or causal mechanisms.
- **Not yet runtime-tested or deployed**; this is a draft for review.
- The package's plot method may need adjustment after live testing.
