"""Cutout preparation, masking and 50-member ensemble inference for the LoTSS S_RG score."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path

import numpy as np
import onnxruntime as ort
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.path import Path as PolygonPath
from PIL import Image
from scipy import ndimage
from skimage import measure, morphology
from skimage.transform import resize

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "ensemble_manifest.json"
MODEL_DIR = Path(os.environ.get("LOFAR_SRG_MODELS", ROOT / "models"))

CLASSES = (
    "Amorphous", "Asymmetric core-jet", "Bright-core", "Centre-bright (FRI)",
    "Core-displaced", "Edge-bright (FRII) linear", "Hybrid (FRX) unclear",
    "Large FRII", "Spiral", "Structured-bent", "Unresolved",
)
RG_CLASSES = (1, 2, 3, 4, 5, 6, 7, 9)
SURVEYS = {"dr2": "https://lofar-surveys.org/dr2-cutout.fits",
           "dr3": "https://lofar-surveys.org/dr3-cutout.fits"}
F_LAS = 1.5
F_CLIP = 3.0
MIN_LAS = 15.0
IMAGE_SIZE = 128
EPS = 1e-10


def fetch_cutout(*, ra: float, dec: float, las: float, survey: str = "dr2") -> fits.PrimaryHDU:
    if survey not in SURVEYS:
        raise ValueError(f"survey must be one of {sorted(SURVEYS)}")
    if not (0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0):
        raise ValueError("RA must lie in [0, 360) and Dec in [-90, 90] degrees")
    if las <= 0:
        raise ValueError("LAS must be positive (arcsec)")
    size_arcmin = max(2.0, F_LAS * max(las, MIN_LAS) / 60.0 + 0.2)
    response = requests.get(SURVEYS[survey], params={"pos": f"{ra},{dec}", "size": f"{size_arcmin:.6f}"},
                            timeout=180)
    if response.status_code == 404:
        raise ValueError(f"Position not covered by LoTSS-{survey.upper()}")
    response.raise_for_status()
    if response.headers.get("content-type") != "application/fits":
        raise ValueError(f"Survey server did not return FITS: {response.text[:200]}")
    with fits.open(io.BytesIO(response.content)) as hdul:
        return fits.PrimaryHDU(data=hdul[0].data.copy(), header=hdul[0].header.copy())


def read_fits(path: str | Path) -> fits.PrimaryHDU:
    with fits.open(path) as hdul:
        hdu = next((h for h in hdul if h.data is not None and h.data.ndim >= 2), None)
        if hdu is None:
            raise ValueError(f"No image data in {path}")
        return fits.PrimaryHDU(data=hdu.data.copy(), header=hdu.header.copy())


def crop(*, hdu: fits.PrimaryHDU, ra: float | None = None, dec: float | None = None,
         las: float | None = None) -> np.ndarray:
    """Square cutout of side F_LAS x LAS centred on (ra, dec); NaNs are filled with the minimum."""
    data = np.squeeze(hdu.data).astype(np.float32)
    if data.ndim != 2:
        raise ValueError(f"Expected a 2-D image, got shape {data.shape}")
    wcs = WCS(hdu.header).celestial
    if (ra is None) != (dec is None):
        raise ValueError("Give both RA and Dec, or neither to use the image centre")
    centre = (SkyCoord(ra, dec, unit="deg") if ra is not None
              else ((data.shape[1] - 1) / 2.0, (data.shape[0] - 1) / 2.0))
    if las is None:
        size = min(data.shape)
    else:
        pixel_arcsec = round(float(proj_plane_pixel_scales(wcs)[0]) * 3600.0, 4)
        size = round(F_LAS * max(las, MIN_LAS) / pixel_arcsec)
    image = Cutout2D(data, centre, size, wcs=wcs, mode="partial").data.copy()
    finite = np.isfinite(image)
    if not finite.any():
        raise ValueError("Cutout contains no valid pixels")
    image[~finite] = image[finite].min()
    return image


def polygon_mask(*, shape: tuple[int, int], polygons: list[list[list[float]]]) -> np.ndarray:
    """Pixels inside any polygon given as [[x, y], ...] in cutout pixel coordinates."""
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    points = np.column_stack((xx.ravel(), yy.ravel()))
    mask = np.zeros(points.shape[0], dtype=bool)
    for polygon in polygons:
        if len(polygon) < 3:
            raise ValueError("Each polygon needs at least three vertices")
        mask |= PolygonPath(np.asarray(polygon, dtype=float)).contains_points(points, radius=0.5)
    return mask.reshape(shape)


def compact_source_mask(*, image: np.ndarray, threshold_sigma: float = 5.0, min_area: int = 5,
                        max_area: int = 300, max_eccentricity: float = 0.9,
                        dilation: int = 2) -> np.ndarray:
    """Bright, small, roundish components (unrelated compact sources), dilated by `dilation` pixels."""
    median = float(np.median(image))
    sigma = max(1.4826 * float(np.median(np.abs(image - median))), 1e-12)
    labels = measure.label(ndimage.gaussian_filter(image, 1.0) > median + threshold_sigma * sigma)
    mask = np.zeros(image.shape, dtype=bool)
    for region in measure.regionprops(labels):
        if min_area <= region.area <= max_area and region.eccentricity <= max_eccentricity:
            mask[labels == region.label] = True
    return morphology.binary_dilation(mask, morphology.disk(dilation)) if dilation and mask.any() else mask


def excluded(*, shape: tuple[int, int], keep: np.ndarray | None = None,
             remove: np.ndarray | None = None) -> np.ndarray:
    """Pixels outside a non-empty `keep` region or inside `remove`."""
    replace = np.zeros(shape, dtype=bool)
    for mask in (keep, remove):
        if mask is not None and mask.shape != shape:
            raise ValueError(f"Mask shape {mask.shape} differs from image shape {shape}")
    if keep is not None and keep.any():
        replace |= ~keep
    if remove is not None:
        replace |= remove
    return replace


def preprocess(*, image: np.ndarray, keep: np.ndarray | None = None,
               remove: np.ndarray | None = None, rounding: bool = False) -> np.ndarray:
    """Clip at F_CLIP x sigma, resize to 128 x 128 and min-max normalise to uint8.

    Pixels outside `keep` or inside `remove` are excluded from the sigma estimate and set to the
    clipping floor. Without masks this reproduces the training images of the catalogue, which were
    quantised by truncation; the thesis relic inputs were quantised by rounding (`rounding=True`).
    """
    replace = excluded(shape=image.shape, keep=keep, remove=remove)
    if replace.all():
        raise ValueError("The mask removes every pixel")
    _, _, std = sigma_clipped_stats(image[~replace], sigma=3.0, maxiters=10)
    floor = F_CLIP * float(std)
    clipped = np.clip(image, floor, np.inf)
    clipped[replace] = floor
    resized = resize(clipped, (IMAGE_SIZE, IMAGE_SIZE))
    lo, hi = float(resized.min()), float(resized.max())
    if not hi > lo:
        raise ValueError("Image is constant after clipping; nothing exceeds the clipping floor")
    scaled = 255 * (resized - lo) / (hi - lo)
    return (np.round(scaled) if rounding else scaled).astype(np.uint8)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_models(*, model_dir: Path = MODEL_DIR) -> None:
    manifest = json.loads(MANIFEST.read_text())
    record = manifest.get("zenodo_record")
    if not record:
        raise RuntimeError(f"No zenodo_record in {MANIFEST}; place the ONNX files in {model_dir} manually")
    model_dir.mkdir(parents=True, exist_ok=True)
    for member in manifest["members"]:
        target = model_dir / member["file"]
        if target.exists() and _sha256(target) == member["sha256"]:
            continue
        print(f"Downloading {member['file']}")
        partial = target.with_suffix(".part")
        with requests.get(f"https://zenodo.org/records/{record}/files/{member['file']}?download=1",
                          stream=True, timeout=300) as response:
            response.raise_for_status()
            with partial.open("wb") as stream:
                for chunk in response.iter_content(1 << 20):
                    stream.write(chunk)
        if _sha256(partial) != member["sha256"]:
            partial.unlink()
            raise RuntimeError(f"Checksum mismatch for {member['file']}")
        partial.replace(target)


class Ensemble:
    def __init__(self, *, model_dir: Path = MODEL_DIR, verify: bool = True) -> None:
        members = json.loads(MANIFEST.read_text())["members"]
        paths = [model_dir / m["file"] for m in members]
        missing = [p.name for p in paths if not p.exists()]
        if missing:
            raise FileNotFoundError(f"{len(missing)} ensemble members missing in {model_dir}; "
                                    "run `python lofar_srg.py --download-models`")
        if verify:
            bad = [m["file"] for m, p in zip(members, paths) if _sha256(p) != m["sha256"]]
            if bad:
                raise ValueError(f"Checksum mismatch: {bad}")
        self.sessions = [ort.InferenceSession(str(p), providers=["CPUExecutionProvider"]) for p in paths]

    def predict(self, *, images: np.ndarray) -> list[dict]:
        """images: uint8 array (N, 128, 128) from `preprocess`."""
        if images.ndim != 3 or images.shape[1:] != (IMAGE_SIZE, IMAGE_SIZE) or images.dtype != np.uint8:
            raise ValueError("Expected uint8 images of shape (N, 128, 128)")
        tensors = (images.astype(np.float32) / 255.0)[:, None]
        logits = np.stack([np.atleast_2d(s.run(None, {"input": tensors})[0]) for s in self.sessions])
        logits = logits.astype(np.float64)
        probs = np.exp(logits - logits.max(axis=2, keepdims=True))
        probs /= probs.sum(axis=2, keepdims=True)
        return [_score(member_probs=probs[:, i]) for i in range(images.shape[0])]


def _score(*, member_probs: np.ndarray) -> dict:
    """S_RG = P_RG * A * (1 - H/log K) * (1 - MI/log K) for member probabilities (M, K)."""
    log_k = np.log(len(CLASSES))
    mean = member_probs.mean(axis=0)
    label = int(mean.argmax())
    entropy = float(-np.sum(mean * np.log(np.clip(mean, EPS, None))))
    member_entropy = float(-np.sum(member_probs * np.log(np.clip(member_probs, EPS, None)), axis=1).mean())
    mutual_info = max(0.0, entropy - member_entropy)
    agreement = float(np.mean(member_probs.argmax(axis=1) == label))
    p_rg = float(mean[list(RG_CLASSES)].sum())
    reliability = agreement * float(np.clip(1 - entropy / log_k, 0, 1)) * float(np.clip(1 - mutual_info / log_k, 0, 1))
    s_rg = p_rg * reliability
    return {"label": label, "class": CLASSES[label], "confidence": float(mean[label]),
            "probabilities": {name: float(p) for name, p in zip(CLASSES, mean)},
            "entropy": entropy, "mutual_info": mutual_info, "agreement": agreement,
            "P_RG": p_rg, "R": reliability, "S_RG": s_rg,
            "band": "high" if s_rg >= 0.75 else "transitional" if s_rg >= 0.45 else "low"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-models", action="store_true", help="Fetch and verify the ONNX ensemble")
    parser.add_argument("--fits", type=Path, help="Local FITS image instead of a LoTSS download")
    parser.add_argument("--ra", type=float, help="Centre RA (deg)")
    parser.add_argument("--dec", type=float, help="Centre Dec (deg)")
    parser.add_argument("--las", type=float, help="Largest angular size (arcsec); cutout side is 1.5 x LAS")
    parser.add_argument("--survey", choices=sorted(SURVEYS), default="dr2")
    parser.add_argument("--polygons", type=Path, help="JSON list of polygons [[x, y], ...] in cutout pixels")
    parser.add_argument("--mode", choices=["keep", "remove"], default="keep",
                        help="keep: polygons enclose the source; remove: polygons cover contaminants")
    parser.add_argument("--auto-mask", action="store_true", help="Remove compact background sources")
    parser.add_argument("--round", action="store_true",
                        help="Round instead of truncate to uint8, as for the thesis relic inputs")
    parser.add_argument("--save-png", type=Path, help="Write the 128 x 128 model input")
    args = parser.parse_args()

    if args.download_models:
        download_models()
        return
    if args.fits is None and None in (args.ra, args.dec, args.las):
        parser.error("give --ra, --dec and --las, or --fits")
    hdu = read_fits(args.fits) if args.fits else fetch_cutout(ra=args.ra, dec=args.dec, las=args.las,
                                                              survey=args.survey)
    image = crop(hdu=hdu, ra=args.ra, dec=args.dec, las=args.las)
    keep = remove = None
    if args.polygons:
        mask = polygon_mask(shape=image.shape, polygons=json.loads(args.polygons.read_text()))
        keep, remove = (mask, None) if args.mode == "keep" else (None, mask)
    if args.auto_mask:
        compact = compact_source_mask(image=image)
        remove = compact if remove is None else remove | compact
    model_input = preprocess(image=image, keep=keep, remove=remove, rounding=args.round)
    if args.save_png:
        Image.fromarray(model_input, mode="L").save(args.save_png)
    print(json.dumps(Ensemble().predict(images=model_input[None])[0], indent=2))


if __name__ == "__main__":
    main()
