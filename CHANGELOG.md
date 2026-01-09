
# 1.6

* `aln_to_base_and_position_table()` now uses data.table library and indexing for a 100x speed boost.

* Updated test for `aln_to_base_and_position_table()`

* `parallel_aln_to_base_and_position_tables()` and corresponding test and documentation were removed. 

* Edit site base editing calculations were turned into a formal function called `characterize_edit_site`

* BEiGUIDE can now processes edit_site objects in parallel. 

* Corrections to function calls in README.

* More usage info added to README.

* Added CHANGELOG.md

