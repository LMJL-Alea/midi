## Resubmission

This is a resubmission. In this version I have:

* Expanded the `Description` field in the `DESCRIPTION` file to better describe 
the package and its purpose;
* Added references both to the `DESCRIPTION` file and the documentation of the
classes;
* Replaced the use of `dontrun{}` in the examples with `donttest{}` to avoid
running the examples during the CRAN checks;
* Added a `seed` option where necessary to avoid fixing the seed in the code.

## Test environments
* local macOS R installation, R 4.3.3
* continuous integration via GH actions:
  * macOS latest release
  * windows latest release
  * ubuntu 20.04 latest release and devel
* [win-builder](https://win-builder.r-project.org/) (release, devel and old-release)
* [macbuilder](https://mac.r-project.org/macbuilder/)
* [R-hub](https://builder.r-hub.io)
  - Windows Server 2022, R-devel, 64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
