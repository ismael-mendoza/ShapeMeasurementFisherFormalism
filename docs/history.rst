Revision History
================

v0.3
----
- Add tutorial notebook.
- Add HSM size and shape analysis for each galaxy (ignoring overlaps)

v0.2
----
- Download from `here <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending/releases/tag/v0.2>`_.
- Refactor code and add sphinx documentation.
- Rename `galsimcat.py` as `simulate.py` and `lsst2wl.py` as `dbquery.py`.
- Add new programs `display.py`, `fisher.py`.
- Fix catalog half-light radius bug: use sqrt(a*b) instead of catalog HLR value.
- Include AGN component when present in the catalog.
- Use Airy optical model instead of increasing atmospheric FWHM by 0.4 arcsec.
- Implement new blending figures of merit.
- Version used for DESC highlight project data products and paper draft v0.2

v0.1
----
- Download from `here <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending/releases/tag/v0.1>`_.
- Version used to prepare results for the `Dec 2013 LSST DESC meeting <https://indico.bnl.gov/conferenceDisplay.py?confId=691>`_ in Pittsburgh.
