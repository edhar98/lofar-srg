# LoTSS radio-galaxy similarity score (S_RG)

Bachelor's thesis *Classification and Similarity Scoring of Extended LOFAR Radio Sources via
Self-Supervised Learning* (Edgar Harutyunyan, Universität Leipzig, September 2026), the resulting
morphology catalogue of 32 514 LoTSS-DR2 sources, and a tool that classifies your own cutouts with
the same 50-member ensemble.

| Path | Content |
|---|---|
| `thesis.pdf` | Thesis as submitted |
| `thesis/` | LaTeX source (`make` builds `main.pdf`; needs pdflatex + bibtex) |
| `catalog_srg.csv` | Catalogue: positions, LAS, ensemble class, class probabilities, uncertainties, S_RG |
| `CLASSES.md` | The 11 morphology classes |
| `lofar_srg.py` | Cutout download, preprocessing, masking, ensemble inference (library + CLI) |
| `app.py` | Browser interface for the same steps |
| `ensemble_manifest.json` | File names, SHA-256 checksums and Zenodo record of the 50 ONNX models |

## Install

Python 3.12:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python lofar_srg.py --download-models    # 50 ONNX files, 2.2 GB, checksum-verified, into models/
```

Set `LOFAR_SRG_MODELS=/path/to/dir` to keep the models elsewhere.

## Browser interface

```bash
python app.py                            # opens http://127.0.0.1:7860
```

1. Enter RA, Dec (deg) and the largest angular size LAS (arcsec) and pick LoTSS-DR2 or DR3,
   or upload your own FITS image. The cutout has side 1.5 x LAS, as for the training images.
2. Optionally mask: paint over contaminating sources (*remove*) or over the source to keep
   (*keep*; everything else is excluded), and/or let the tool remove compact background sources.
3. **Run model** shows the excluded pixels, the 128 x 128 model input, the class probabilities
   and S_RG with its components.

## Command line

```bash
python lofar_srg.py --ra 15.974684 --dec 19.265668 --las 170.2                    # LoTSS-DR2
python lofar_srg.py --ra 90.8333 --dec 42.3083 --las 474 --survey dr3 --auto-mask
python lofar_srg.py --fits my_image.fits --ra .. --dec .. --las .. \
                    --polygons polygons.json --mode keep --save-png input.png
```

`polygons.json` is a list of polygons `[[x, y], ...]` in pixel coordinates of the 1.5 x LAS cutout.
Output is JSON with the class probabilities, entropy, mutual information, agreement, P_RG, R and S_RG.

## Method

Preprocessing is identical to the training images: square cutout of side 1.5 x max(LAS, 15")
around the given position, clipping at 3 sigma (sigma-clipped standard deviation), resizing to
128 x 128 px and min-max normalisation to 8 bit. Masked pixels are excluded from the sigma estimate
and set to the clipping floor, as in the relic study of the thesis.

For member probabilities p_m over K = 11 classes with mean p, predicted class c = argmax p:

```text
P_RG  = sum of p over the radio-galaxy classes {1,2,3,4,5,6,7,9}
H     = entropy of p;  MI = H - mean member entropy;  A = fraction of members predicting c
R     = A (1 - H/log K) (1 - MI/log K)
S_RG  = P_RG R        high >= 0.75, transitional 0.45-0.75, low < 0.45
```

S_RG measures morphological similarity to the radio-galaxy classes learned from LoTSS-DR2; it is
not a probability of physical source type. See the thesis for training, validation and limitations.

## Reproducibility

- Catalogue sources downloaded through this tool give the training images to within one 8-bit
  level and the catalogue labels and S_RG to within 1e-4.
- The catalogue images of sources with LAS > 240" (5 % of the sample) were made from 6' downloads,
  so their outer parts were padded; this tool always downloads the full 1.5 x LAS footprint and
  can therefore score these sources differently from the catalogue.
- The relic inputs in the thesis were quantised to 8 bit by rounding instead of truncation. Use
  `--round` to reproduce the relic table; the score is sensitive to this at the 0.01-0.06 level.

## Data and credits

Images are from the LOFAR Two-metre Sky Survey (LoTSS) DR2 (Shimwell et al. 2022, A&A 659, A1) and
DR3 (Shimwell et al. 2026, A&A 707, A198), retrieved from <https://lofar-surveys.org>; LAS values are from the
DR2 value-added catalogue (Hardcastle et al. 2023, A&A 678, A151). The encoder, contrastive
framework and classifier pipeline are adapted from
[RadioSimCLR](https://github.com/NicolasBaronPerez/RadioSimCLR) and
[ClassifierFR](https://github.com/NicolasBaronPerez/ClassifierFR) by N. Baron Perez.

Code is released under the MIT licence (`LICENSE`). Third-party figures in the thesis are credited
in their captions and remain under their original terms.
