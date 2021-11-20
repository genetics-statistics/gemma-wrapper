
## ChangeLog v0.99.4 (2021xxxx)

* Fixed problem of incomplete cache files
* Adding locking so the exact same command does not run twice - with a risk of overwriting output files

## ChangeLog v0.99.3 (20210822)

* Fetch chromosome names from the annotation file - working of --loco switch changed!
* Use --chromosomes switch to override which chromosomes to compute
* Also allow for parallel GWA compute of non-LOCO (unLOCO?)
* Added metrics page

## ChangeLog v0.99.2 (20210808)

* Use isolated tmpdir for GEMMA and don't overwrite output files if
  they appeared during run (make it more transactional)
* Run GEMMA in parallel for LOCO with 5x speedups
