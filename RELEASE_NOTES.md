## ChangeLog v0.99.7 ()

* Less aggressive trials of parallel compute - should really scale to size
* Check for JSON switch

## ChangeLog v0.99.6 (20220122)

Important bug fix release!!

* Fixed bug that had creeped in and always returned the same hash for all runs (terrible)
* Zoomed in on tests to make sure this won't happen again
* Found another potential disaster with too long a command line for parallel jobs

## ChangeLog v0.99.5

* Move parallel joblog out of the way after a run

## ChangeLog v0.99.4 (20211125)

* Fixed problem of incomplete cache files by introducing `parallel --joblog`
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
