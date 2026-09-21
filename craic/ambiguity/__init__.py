"""Tools for finding, scoring, and resolving ambiguously aligned regions.

Four complementary views of alignment uncertainty:
    reliability   - per-column / per-residue confidence + live masking
    disagreement  - where independent aligners disagree
    sandbox       - human-in-the-loop local realignment
    posterior     - the pair-HMM posterior itself, as alternative homologies
"""
