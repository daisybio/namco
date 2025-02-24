# v0.2.2 (2023-10-08)

* Made sure that Phyloseq objects with a single sample are displayed without an error.
* Fixed snapshot path handling (broken in some cases)

# v0.2.1 (2022-11-21)

* Fixed snapshots on windows (images were empty)
* `krona_opts` are applied in interactive charts as well with the default settings

# v0.2 (2022-11-15)

This release brings a whole lot of improvements.

**Main improvements**

* The restriction to phyloseq as input data is removed, there
  are now several ways of providing input data
* Custom coloring improved, more options
* Temporary files are cleaned up in most cases (except if used from the console)
* The special Perl import script is not needed anymore
* A tutorial was added
* Rendering now works on Windows
 
**Other changes** (from commit log)

* Captions are now possible in HTML output (with bookdown)
* Abundance summary functions can be configured
* snapshot_format is auto-recognized based  on output file extension
* Krona Javascript can be included from remote sources and
  is minified by default to save space

# v0.1 (2022-09-20)

Initial release
