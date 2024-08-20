
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# Peripheral immune cells contributing to neurodegeneration

## Background

<!-- Add a description of your project. -->

This project explores whether peripheral immune cells mediate
neurodegenerative diseases in an ancestry specific manner. It does this
in 2 ways:

- First it explores the extent and pattern of expression of
  neurodegenerative diseases in peripheral immune cells
- It then uses genetic colocalisation to explore whether the expression
  of these genes in peripheral immune cells is causally linked with the
  risk of disease.

Noting that regional selection pressures will affect environmental-host
interactions in a region specific manner, this analysis was undertaken
comparing African, East Asian and European ancestries.

## Code contents

Within this repository you will find:

| Directory                        | Description                                                                                                                                                |
|----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [docs](docs)                     | Contains all `.Rmd`s and their corresponding `.html`s describing analyses performed for this project. These can be view interactively at: [link](#TODO)    |
| [logs](logs)                     | For any scripts that were run outside of an `.Rmd` (e.g. scripts from the [scripts](scripts) directory), a log file was recorded and can be accessed here. |
| [raw_data](raw_data)             | External tables and genomic datasets used in analyses (not linked to github due to size, see published article for sources).                               |
| [processed_data](processed_data) | Results output from colocalisation analysis pipelines (not linked to github due to size).                                                                  |
| [derived_data](derived_data)     | Processed results.                                                                                                                                         |
| [renv](renv)                     | `renv`-related scripts                                                                                                                                     |
| [r_scripts](r_scripts)           | Contains R based scripts used throughout the analysis, also referenced in respective `.Rmd`.                                                               |
| [scripts](scripts)               | Contains analysis scripts.                                                                                                                                 |

## Reproducibility

<!-- Modify selection below depending on how package dependencies have been managed. -->

### `renv`

<!-- Consider using renv for reproducibility. Delete this section if you will not be doing this. -->

This repository uses [`renv`](https://rstudio.github.io/renv/index.html)
to create a reproducible environment for this R project.

1.  When you first launches this project, `renv` should automatically
    bootstrap itself, thereby downloading and installing the appropriate
    version of `renv` into the project library.
2.  After this has completed, you can use `renv::restore()` to restore
    the project library locally on your machine.

For more information on collaborating with `renv`, please refer to this
[link](https://rstudio.github.io/renv/articles/collaborating.html).

## License

<!-- For analyses, an MIT license can be added to the project using usethis::use_mit_license(). -->
<!-- If you don't end up using an MIT license, edit below. -->

The code in this repository is released under an MIT license. This
repository is distributed in the hope that it will be useful to the
wider community, but without any warranty of any kind. Please see the
[LICENSE](LICENSE) file for more details.

## Citation

<!-- Add any necessary software citations -->

This package makes use of the [coloc package, version
5.1.0.1](https://cran.r-project.org/package=coloc), as well as wrapper
functions found in [colochelper version
0.99.1](http://dx.doi.org/10.5281/zenodo.5011869).
