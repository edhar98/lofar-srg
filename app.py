"""Gradio interface: make a LoTSS cutout, mask it, and score it with the 50-member ensemble."""

from __future__ import annotations

import gradio as gr
import numpy as np
from skimage.transform import resize

import lofar_srg as L

CANVAS = 512
SOURCE_LOTSS = "LoTSS download"
SOURCE_UPLOAD = "Upload FITS"
MODE_KEEP = "Keep painted region (paint over the source)"
MODE_REMOVE = "Remove painted regions (paint over contaminants)"
_ensemble: L.Ensemble | None = None


def get_ensemble() -> L.Ensemble:
    global _ensemble
    if _ensemble is None:
        try:
            _ensemble = L.Ensemble()
        except FileNotFoundError:
            L.download_models()
            _ensemble = L.Ensemble()
    return _ensemble


def display(image: np.ndarray, overlay: np.ndarray | None = None) -> np.ndarray:
    """North-up RGB rendering with a 1-99.7 percentile stretch; `overlay` pixels are tinted red."""
    lo, hi = np.percentile(image, [1.0, 99.7])
    grey = np.clip((image - lo) / (hi - lo if hi > lo else 1.0), 0.0, 1.0)
    rgb = np.repeat(grey[..., None], 3, axis=2)
    if overlay is not None:
        rgb[overlay] = 0.6 * rgb[overlay] + 0.4 * np.array([1.0, 0.0, 0.0])
    rgb = resize(rgb, (CANVAS, CANVAS), order=1)
    return (255 * np.flipud(rgb)).astype(np.uint8)


def painted_mask(editor: dict | None, shape: tuple[int, int]) -> np.ndarray:
    layers = [layer for layer in (editor or {}).get("layers") or [] if layer is not None]
    if not layers:
        return np.zeros(shape, dtype=bool)
    painted = np.any([layer[..., 3] > 0 for layer in layers], axis=0)
    return np.flipud(resize(painted.astype(float), shape, order=0, anti_aliasing=False) > 0.5)


def make_cutout(source: str, ra: float | None, dec: float | None, las: float | None, survey: str,
                upload: str | None) -> tuple[dict, np.ndarray, str]:
    try:
        if source == SOURCE_LOTSS:
            if None in (ra, dec, las):
                raise ValueError("RA, Dec and LAS are required for a LoTSS download")
            hdu = L.fetch_cutout(ra=ra, dec=dec, las=las, survey=survey)
        else:
            if upload is None:
                raise ValueError("Upload a FITS image first")
            hdu = L.read_fits(upload)
            ra, dec = (None, None) if not ra and not dec else (ra, dec)
            las = las or None
        image = L.crop(hdu=hdu, ra=ra, dec=dec, las=las)
    except Exception as exc:
        raise gr.Error(str(exc)) from exc
    view = display(image)
    status = f"Cutout: {image.shape[1]} x {image.shape[0]} px. Paint on the image to mask, then press **Run model**."
    return {"background": view, "layers": [], "composite": view}, image, status


def run_model(editor: dict | None, image: np.ndarray | None, mode: str, auto: bool, threshold: float,
              max_area: float) -> tuple[np.ndarray, np.ndarray, dict, str]:
    if image is None:
        raise gr.Error("Make a cutout first")
    painted = painted_mask(editor, image.shape)
    keep, remove = (painted, None) if mode == MODE_KEEP else (None, painted)
    if auto:
        compact = L.compact_source_mask(image=image, threshold_sigma=threshold, max_area=int(max_area))
        remove = compact if remove is None else remove | compact
    try:
        model_input = L.preprocess(image=image, keep=keep, remove=remove)
        result = get_ensemble().predict(images=model_input[None])[0]
    except Exception as exc:
        raise gr.Error(str(exc)) from exc
    replaced = L.excluded(shape=image.shape, keep=keep, remove=remove)
    metrics = (
        f"### S_RG = {result['S_RG']:.3f} ({result['band']})\n\n"
        f"Predicted class: **{result['class']}** ({result['confidence']:.1%})\n\n"
        "| P_RG | R | agreement | entropy | mutual information |\n|---|---|---|---|---|\n"
        f"| {result['P_RG']:.3f} | {result['R']:.3f} | {result['agreement']:.0%} | "
        f"{result['entropy']:.3f} | {result['mutual_info']:.3f} |\n\n"
        "Bands: high >= 0.75, transitional 0.45-0.75, low < 0.45. "
        "S_RG is a morphological similarity index, not a probability of physical source type."
    )
    return display(image, overlay=replaced), np.flipud(model_input), result["probabilities"], metrics


with gr.Blocks(title="LoTSS radio-galaxy similarity") as demo:
    gr.Markdown(
        "# LoTSS radio-galaxy similarity (S_RG)\n"
        "Make a LoTSS cutout of side 1.5 x LAS, optionally mask it, and classify it with the "
        "50-member ensemble into 11 morphology classes. Masked pixels are excluded from the noise "
        "estimate and set to the 3-sigma clipping floor before resizing to 128 x 128 px."
    )
    crop_state = gr.State(None)
    with gr.Row():
        with gr.Column(scale=1):
            source = gr.Radio([SOURCE_LOTSS, SOURCE_UPLOAD], value=SOURCE_LOTSS, label="Input")
            ra = gr.Number(label="RA (deg)", value=None)
            dec = gr.Number(label="Dec (deg)", value=None)
            las = gr.Number(label="Largest angular size, LAS (arcsec)", value=None)
            survey = gr.Radio(sorted(L.SURVEYS), value="dr2", label="Survey (download only)")
            upload = gr.File(label="FITS image (upload only; RA = Dec = 0 uses the image centre, "
                                   "LAS = 0 the whole image)", type="filepath")
            cut_button = gr.Button("Make cutout", variant="primary")
            status = gr.Markdown()
        with gr.Column(scale=2):
            editor = gr.ImageEditor(type="numpy", sources=(), transforms=(), layers=False,
                                    brush=gr.Brush(colors=["#ff3030"], color_mode="fixed"),
                                    canvas_size=(CANVAS, CANVAS), height=620,
                                    label="Cutout (north up): paint to mask")
            mode = gr.Radio([MODE_KEEP, MODE_REMOVE], value=MODE_REMOVE, label="Painted pixels")
            with gr.Row():
                auto = gr.Checkbox(label="Also remove compact background sources", value=False)
                threshold = gr.Slider(3, 15, value=5, step=0.5, label="Detection threshold (sigma)")
                max_area = gr.Slider(20, 2000, value=300, step=10, label="Max compact-source area (px)")
            run_button = gr.Button("Run model", variant="primary")
    with gr.Row():
        mask_view = gr.Image(label="Excluded pixels (red)", height=320)
        input_view = gr.Image(label="Model input (128 x 128)", height=320, image_mode="L")
        probabilities = gr.Label(label="Class probabilities", num_top_classes=5)
    metrics = gr.Markdown()
    gr.Examples(
        examples=[[SOURCE_LOTSS, 15.974684, 19.265668, 170.2, "dr2"],
                  [SOURCE_LOTSS, 90.8333333, 42.3083333, 474.0, "dr3"]],
        inputs=[source, ra, dec, las, survey],
        label="Examples: a catalogue Large FRII (DR2) and the Toothbrush relic B (DR3)",
    )
    cut_button.click(make_cutout, inputs=[source, ra, dec, las, survey, upload],
                     outputs=[editor, crop_state, status])
    run_button.click(run_model, inputs=[editor, crop_state, mode, auto, threshold, max_area],
                     outputs=[mask_view, input_view, probabilities, metrics])

if __name__ == "__main__":
    demo.launch()
