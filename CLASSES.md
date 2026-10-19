# Morphology Classes

This document describes the 11 radio galaxy morphology classes used by the classifier.
Class labels map to the `ensemble_label` / `ensemble_class` columns in `catalog_srg.csv`.

| ID | Name | Count in catalog |
|----|------|-----------------|
| 0  | Amorphous | 390 |
| 1  | Asymmetric core-jet (ACJ) | 1 976 |
| 2  | Bright-core | 3 728 |
| 3  | Centre-bright (FRI) | 3 323 |
| 4  | Core-Displaced | 980 |
| 5  | Edge-bright (FRII) linear | 7 574 |
| 6  | Hybrid (FRX) unclear | 2 865 |
| 7  | Large FRII | 2 425 |
| 8  | Spirals | 963 |
| 9  | Structured-Bent | 7 197 |
| 10 | Unresolved | 1 093 |

---

## Class descriptions

### 0 — Amorphous
Diffuse, extended emission with no clearly defined jets, lobes, or hotspots. The radio structure
lacks the organised geometry expected of active jets. These sources often represent remnant or
restarted AGN where the plasma has expanded and faded after jet activity ceased, or genuinely
relaxed sources whose morphology has been disrupted by the environment.

### 1 — Asymmetric core-jet (ACJ)
A bright central core with a single visible jet extending in one direction. The counter-jet is
below the detection threshold, either due to Doppler boosting (the approaching jet is amplified,
the receding jet suppressed) or asymmetric absorption. Common in AGN whose jets are oriented at
a moderate angle to the line of sight, or in genuinely one-sided sources in dense environments.

### 2 — Bright-core
Core-dominated sources where the unresolved or barely resolved nucleus contributes the majority
of the total flux. Extended emission (jets or lobes) is faint or absent at LoTSS resolution.
Includes compact steep-spectrum and gigahertz-peaked sources, as well as sources where the
extended emission has faded while the core remains active.

### 3 — Centre-bright (FRI)
Fanaroff–Riley Type I morphology. Twin jets emerge from the nucleus, are brightest close to the
core, and fade outward into diffuse lobes without a compact terminal hotspot (edge-darkened).
FRI sources are typically found in lower radio-power AGN hosted by early-type galaxies in group
or cluster environments, where jet deceleration by the denser interstellar medium produces the
characteristic brightness gradient.

### 4 — Core-Displaced
Sources where the radio core does not coincide with the geometric centroid of the emission.
One lobe or jet dominates strongly over the other, or the core is offset from the midpoint of
the source. Possible causes include strongly asymmetric environments (one jet propagates into
a denser medium), jet precession, or past merger activity that has displaced the host nucleus
relative to older radio plasma. A fraction of these sources may also be intrinsically
core-dominated (Bright-core, class 2) where the catalogued RA/DEC position in LoTSS is placed
slightly off the true peak, causing the bright nucleus to appear displaced from the image centre
in the 128 × 128 cutout.

### 5 — Edge-bright (FRII) linear
Fanaroff–Riley Type II morphology. Collimated twin jets terminate in compact, high-surface-brightness
hotspots at the outer edges of the lobes (edge-brightened). Backflowing plasma fills the lobes
behind the hotspots. FRII sources span a wide range of angular sizes and are typically found in
higher radio-power AGN. The linear qualifier indicates that the two lobes and core are well-aligned
along a single axis.

### 6 — Hybrid (FRX) unclear
Sources that display mixed or ambiguous morphology — for example, FRI-like emission on one side
and FRII-like on the other (HyMoRS, Hybrid Morphology Radio Sources), or sources where the
available resolution and sensitivity make a definitive FRI/FRII classification uncertain. The FRX
label is used here to denote this transitional or unclear category.

### 7 — Large FRII
FRII morphology at unusually large angular extent (typically LAS ≫ 200 arcsec in this survey).
These are candidate giant radio galaxies (GRGs, projected physical size ≥ 1 Mpc). The hotspot
and lobe structure is consistent with classical FRII, but the source spans a substantial fraction
of the cutout or extends near its edges. Separated from the standard FRII class because the very
large size affects the appearance in the 128 × 128 pixel cutout.

### 8 — Spirals
Radio emission tracing the disk and spiral arms of a star-forming galaxy rather than an AGN jet.
The synchrotron radiation originates from relativistic electrons accelerated by supernova remnants
distributed across the galactic disk. The morphology follows the optical spiral structure and
can superficially resemble bent AGN jets but lacks a bright core or hotspots.

### 9 — Structured-Bent
Sources with clearly bent or curved jet and lobe morphology. Includes Wide-Angle Tailed (WAT)
galaxies, Narrow-Angle Tailed (NAT) / head-tail sources, and C-shaped or S-shaped sources.
Bending arises primarily from ram pressure as the host galaxy moves through the intracluster
medium, making this class a reliable indicator of cluster membership. S-shaped morphologies may
additionally reflect jet precession or binary SMBH activity.

### 10 — Unresolved
Point-like or marginally resolved sources at the LoTSS resolution (∼6 arcsec FWHM). Despite
having a catalog LAS ≥ 60 arcsec (the selection threshold), these sources appear unresolved in
the 128 × 128 pixel cutout, which may indicate an overestimated LAS in the source catalogue,
a compact steep-spectrum morphology, or a source at high redshift where the physical extent
falls below the angular resolution limit.
