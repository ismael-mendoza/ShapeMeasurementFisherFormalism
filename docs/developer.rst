Information for Developers
==========================

Build the documentation
-----------------------

To build a local copy of the HTML documentation, use::

	cd .../descwl/docs
	make html

To view the most recent build of the HTML documentation, point your browser to `.../docs/_build/html/index.html`

To create a tarball snapshot `.../descwl/docs/_build/descwl.tgz` that can be installed on a web server, use::

	cd .../descwl/docs/_build/html
	tar -zcf ../descwl.tgz .

Add a new package module
------------------------

Create a new file `analysis/xxx.py` for module `analysis.xxx` with an initial descriptive docstring.

Add the line `import xxx` to `analysis/__init__.py`.

Add the line `analysis.xxx` to the list of submodules in `docs/src/analysis.rst`.

Create a new documentation file `docs/src/analysis.xxx.rst` containing (replace `xxx` in two places)::

	analysis.xxx module
	=================

	.. automodule:: analysis.xxx
	    :members:
	    :undoc-members:
	    :show-inheritance:

Update the version
------------------

Update the `version` and `release` values in the sphinx configuration file `docs/conf.py`, then commit, tag, and push the new version, e.g.::

	git tag -l v0.2

