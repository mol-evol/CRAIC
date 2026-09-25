"""Color schemes for residues, reliability scores, and posteriors."""

from __future__ import annotations

from PySide6.QtGui import QColor

from ..domain import Alphabet, Level, translate_codon

# Background panel palette (dark, easy on the eyes for long sessions)
BG = QColor("#14181d")
PANEL = QColor("#1b2027")
GRID = QColor("#2a313a")
TEXT = QColor("#dfe6ee")
MUTED = QColor("#8a97a6")
GAP = QColor("#222831")
MASKED = QColor("#39414b")      # a residue the user has masked: kept, but not exported

_NT = {
    "A": "#4CAF50", "C": "#2196F3", "G": "#FF9800", "T": "#E53935",
    "U": "#E53935", "N": "#6b7682", "-": "#222831", ".": "#222831",
}

_DEFAULT = "#3a424c"


def _group_map(groups):
    d = {}
    for residues, hexc in groups:
        for c in residues:
            d[c] = hexc
    return d


# Clustal-like physicochemistry (default)
_CLUSTAL = _group_map([
    ("AVLIMFWY", "#2d7dd2"), ("G", "#f4a259"), ("P", "#e8c547"),
    ("STNQ", "#3bb273"), ("KRH", "#e15554"), ("DE", "#9d4edd"), ("C", "#e0b1cb"),
])
# Zappo (physicochemical)
_ZAPPO = _group_map([
    ("ILVAM", "#ff8c8c"), ("FWY", "#ffb300"), ("KRH", "#5b6bff"),
    ("DE", "#ff2b2b"), ("STNQ", "#26d07c"), ("PG", "#d36bd3"), ("C", "#ffe000"),
])
# Taylor (1997)
_TAYLOR = {
    "A": "#ccff00", "R": "#0000ff", "N": "#cc00ff", "D": "#ff0000", "C": "#ffff00",
    "Q": "#ff00cc", "E": "#ff0066", "G": "#ff9900", "H": "#0066ff", "I": "#66ff00",
    "L": "#33ff00", "K": "#6600ff", "M": "#00ff00", "F": "#00ff66", "P": "#ffcc00",
    "S": "#ff3300", "T": "#ff6600", "W": "#00ccff", "Y": "#00ffcc", "V": "#99ff00",
}
# Hydrophobicity (Kyte-Doolittle): red = hydrophobic, blue = hydrophilic
_KD = {"I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
       "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
       "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5}


def _hydro_hex(v: float) -> str:
    t = (v + 4.5) / 9.0
    if t < 0.5:
        f = t / 0.5
        r = g = int(255 * f); b = 255
    else:
        f = (t - 0.5) / 0.5
        r = 255; g = b = int(255 * (1 - f))
    return f"#{r:02x}{g:02x}{b:02x}"


_HYDRO = {a: _hydro_hex(v) for a, v in _KD.items()}

_AA_SCHEMES = {"clustal": _CLUSTAL, "zappo": _ZAPPO, "taylor": _TAYLOR, "hydrophobicity": _HYDRO}
for _sch in _AA_SCHEMES.values():
    _sch.setdefault("-", "#222831")
    _sch.setdefault("X", "#5a6470")
    _sch.setdefault("*", "#111418")

#: The selectable residue-colouring schemes. Note that all four are *amino-acid*
#: schemes — see :func:`scheme_applies`.
SCHEME_LABELS = ["Clustal", "Zappo", "Taylor", "Hydrophobicity"]


def scheme_applies(level: Level, alphabet: Alphabet) -> bool:
    """Does the chosen scheme affect what is currently on screen?

    Every scheme in :data:`SCHEME_LABELS` colours amino acids; nucleotides are
    coloured by base and ignore the scheme entirely (see :func:`residue_color`,
    whose branches this mirrors). So in a plain nucleotide view the control is
    inert, and offering it there invites the user to change a setting that does
    nothing.
    """
    if alphabet == Alphabet.PROTEIN:
        return True
    return level in (Level.CODON, Level.AA)


def _nt_color(ch: str) -> QColor:
    return QColor(_NT.get(ch.upper(), _DEFAULT))


def _aa_color(ch: str, scheme: str = "clustal") -> QColor:
    return QColor(_AA_SCHEMES.get(scheme, _CLUSTAL).get(ch.upper(), _DEFAULT))


def residue_color(cell: str, level: Level, alphabet: Alphabet,
                  table: int = 1, scheme: str = "clustal") -> QColor:
    """Colour a display cell by the viewed level, under the chosen scheme.

    * nucleotide view  -> colour by base
    * codon view       -> colour the codon by the amino acid it encodes
    * amino-acid view  -> colour by residue
    Amino-acid colours follow ``scheme`` (clustal / zappo / taylor / hydrophobicity).
    """
    if alphabet == Alphabet.PROTEIN:
        return _aa_color(cell, scheme)
    if level == Level.CODON:
        return _aa_color(translate_codon(cell, table), scheme)
    if level == Level.AA:
        return _aa_color(cell, scheme)
    return _nt_color(cell)


def _lerp(a: int, b: int, t: float) -> int:
    return int(round(a + (b - a) * t))


def score_color(v: float) -> QColor:
    """Reliability/agreement ramp: red (0) -> amber (0.5) -> green (1)."""
    if v != v:  # nan
        return QColor("#3a424c")
    v = max(0.0, min(1.0, float(v)))
    if v < 0.5:
        t = v / 0.5
        r, g, b = 0xE5, _lerp(0x39, 0xC1, t), 0x35
    else:
        t = (v - 0.5) / 0.5
        r, g, b = _lerp(0xE5, 0x4C, t), _lerp(0xC1, 0xAF, t), _lerp(0x35, 0x50, t)
    return QColor(r, g, b)


# Perceptual colormaps for the posterior heatmap. Each is a list of RGB anchor
# stops at evenly spaced positions in [0, 1]; we linearly interpolate between the
# two nearest stops. The first few are the standard scientific maps.
_COLORMAPS = {
    "Viridis":   [(68, 1, 84), (59, 82, 139), (33, 144, 140), (93, 201, 99), (253, 231, 37)],
    # Magma trends pink/salmon through the mids; Inferno trends orange/gold — the
    # two only converge at the extremes, so they need enough stops to stay distinct.
    "Magma":     [(0, 0, 4), (28, 16, 68), (79, 18, 123), (129, 37, 129), (181, 54, 122),
                  (229, 80, 100), (251, 135, 97), (254, 194, 135), (252, 253, 191)],
    "Inferno":   [(0, 0, 4), (31, 12, 72), (85, 15, 109), (136, 34, 106), (186, 54, 85),
                  (227, 89, 51), (249, 140, 10), (249, 201, 50), (252, 255, 164)],
    "Turbo":     [(48, 18, 59), (54, 125, 197), (43, 200, 128), (173, 220, 46),
                  (249, 148, 40), (122, 4, 3)],
    "Cividis":   [(0, 32, 76), (45, 87, 125), (124, 123, 120), (190, 164, 99), (255, 233, 69)],
    "Cyan":      [(18, 24, 32), (0, 150, 170), (234, 246, 255)],
    "Greyscale": [(20, 20, 20), (245, 245, 245)],
}
COLORMAP_NAMES = list(_COLORMAPS)


def colormap_color(v: float, name: str = "Viridis") -> QColor:
    """Map ``v`` in [0, 1] to a colour on the named colormap."""
    if v != v:  # nan
        return QColor("#3a424c")
    v = max(0.0, min(1.0, float(v)))
    stops = _COLORMAPS.get(name) or _COLORMAPS["Viridis"]
    n = len(stops) - 1
    pos = v * n
    k = int(pos)
    if k >= n:
        return QColor(*stops[n])
    t = pos - k
    a, b = stops[k], stops[k + 1]
    return QColor(_lerp(a[0], b[0], t), _lerp(a[1], b[1], t), _lerp(a[2], b[2], t))
