"""CRAIC — Conserved-Region Alignment by Iterative Convergence.

A multiple sequence alignment workbench built for ambiguous regions.

The package is intentionally light to import: the GUI and heavy engines are
imported lazily by their own modules so that headless / scripted use never
pays for PySide6.
"""

__version__ = "0.5.9"
