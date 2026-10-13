#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_extract_worldclim_bioclim.py
======================================
Author: Iris Sammarco
Date: 20/07/2026

Purpose
-------
Extract WorldClim v2.1 historical baseline bioclimatic variables (BIO1-BIO19)
for point coordinates stored in a CSV file.

This script is designed to be:
- reproducible,
- transparent,
- cache-aware,
- and easier to audit.

What it does
------------
1. Reads an input CSV with point coordinates.
2. Validates latitude/longitude columns and coordinate ranges.
3. Downloads WorldClim v2.1 bioclim archives only if missing.
4. Extracts raster files into a cache directory only if needed.
5. Samples BIO1-BIO19 for each point using nearest or bilinear interpolation.
6. Writes one output CSV per resolution.
7. Writes a metadata text file per resolution describing what was done.

Important notes
---------------
- The script does NOT add explanatory BIO name columns to the output CSV.
- The script does NOT re-download files if the ZIP already exists.
- By default, it uses nearest-neighbour sampling because that is the most
  reproducible when matching precomputed values from another workflow.
- WorldClim temperature-related variables are commonly scaled; do not assume
  raw values are literal degrees Celsius unless you confirm the stored unit
  convention for the specific raster files you are using.

Input assumptions
-----------------
- CSV is comma-delimited.
- Coordinates are decimal degrees.
- Latitude column name and longitude column name are configurable.
- Input coordinates are WGS84 / EPSG:4326.

Output columns
--------------
Original input columns are preserved.
Added columns:
- bio1 ... bio19

Run with
mamba activate worldclim_env
python 01_extract_worldclim_bioclim.py

Dataset reference
-----------------
WorldClim version 2.1 historical baseline (1970-2000):
Fick, S. E., & Hijmans, R. J. (2017).
WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas.
International Journal of Climatology, 37(12), 4302–4315.
https://doi.org/10.1002/joc.5086
"""

from __future__ import annotations

import csv
import sys
import zipfile
from pathlib import Path

import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
import requests

# ---------------------------------------------------------------------------
# USER CONFIGURATION
# ---------------------------------------------------------------------------

INPUT_CSV = Path(
    "/Anastatica/climate/extract_coordinates/Anastatica_all_coordinates.csv"
)

OUTPUT_DIR = Path(
    "/Anastatica/climate/extract_coordinates"
)

CACHE_DIR = Path(
    "/Anastatica/climate/extract_coordinates/worldclim_cache"
)

LATITUDE_COLUMN = "latitude"
LONGITUDE_COLUMN = "longitude"

# Use "nearest" if you want to match another extraction workflow as closely as
# possible. Use "bilinear" only if you intentionally want interpolated values.
RESAMPLING_METHOD = "nearest"

# Set to True only if you want the script to reuse existing extracted TIFFs.
# If False, it will still avoid re-downloading ZIPs, but will re-check extraction.
REUSE_EXTRACTED_TIFS = True

WORLDCLIM_URLS = {
    "30s": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip",
    "2.5m": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip",
    "5m": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_bio.zip",
    "10m": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip",
}

BIO_CODES = list(range(1, 20))

# ---------------------------------------------------------------------------
# LOGGING HELPERS
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)

# ---------------------------------------------------------------------------
# VALIDATION
# ---------------------------------------------------------------------------

def validate_configuration() -> None:
    if RESAMPLING_METHOD not in {"nearest", "bilinear"}:
        raise ValueError("RESAMPLING_METHOD must be 'nearest' or 'bilinear'.")

def validate_input(df: pd.DataFrame) -> pd.DataFrame:
    if LATITUDE_COLUMN not in df.columns:
        raise KeyError(f"Missing latitude column: {LATITUDE_COLUMN}")
    if LONGITUDE_COLUMN not in df.columns:
        raise KeyError(f"Missing longitude column: {LONGITUDE_COLUMN}")

    out = df.copy()
    out[LATITUDE_COLUMN] = pd.to_numeric(out[LATITUDE_COLUMN], errors="coerce")
    out[LONGITUDE_COLUMN] = pd.to_numeric(out[LONGITUDE_COLUMN], errors="coerce")

    bad = out[out[LATITUDE_COLUMN].isna() | out[LONGITUDE_COLUMN].isna()]
    if not bad.empty:
        raise ValueError(
            f"Found {len(bad)} rows with non-numeric coordinates in "
            f"'{LATITUDE_COLUMN}'/'{LONGITUDE_COLUMN}'."
        )

    bad_range = out[
        (out[LATITUDE_COLUMN] < -90) | (out[LATITUDE_COLUMN] > 90) |
        (out[LONGITUDE_COLUMN] < -180) | (out[LONGITUDE_COLUMN] > 180)
    ]
    if not bad_range.empty:
        raise ValueError(
            f"Found {len(bad_range)} rows with coordinates outside valid WGS84 bounds."
        )

    return out

# ---------------------------------------------------------------------------
# DOWNLOAD AND EXTRACTION
# ---------------------------------------------------------------------------

def download_file(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"Downloading: {url}")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    fh.write(chunk)

def expected_tifs_for_resolution(resolution: str) -> list[Path]:
    extract_dir = CACHE_DIR / f"wc2.1_{resolution}_bio"
    return [extract_dir / f"wc2.1_{resolution}_bio_{bio}.tif" for bio in BIO_CODES]

def ensure_worldclim_data(resolution: str) -> list[Path]:
    url = WORLDCLIM_URLS[resolution]
    zip_name = url.split("/")[-1]
    zip_path = CACHE_DIR / zip_name
    extract_dir = CACHE_DIR / f"wc2.1_{resolution}_bio"

    tif_paths = expected_tifs_for_resolution(resolution)
    if REUSE_EXTRACTED_TIFS and all(p.exists() for p in tif_paths):
        log(f"[{resolution}] Reusing extracted TIFFs from cache.")
        return tif_paths

    if not zip_path.exists():
        download_file(url, zip_path)
    else:
        log(f"[{resolution}] ZIP already exists, skipping download.")

    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zf:
        members = [m for m in zf.namelist() if m.lower().endswith(".tif")]
        if len(members) != 19:
            raise FileNotFoundError(
                f"[{resolution}] Expected 19 TIFFs in archive, found {len(members)}."
            )
        for member in members:
            target = extract_dir / Path(member).name
            if target.exists() and REUSE_EXTRACTED_TIFS:
                continue
            zf.extract(member, extract_dir)
            extracted = extract_dir / member
            if extracted != target:
                target.parent.mkdir(parents=True, exist_ok=True)
                extracted.replace(target)

    tif_paths = expected_tifs_for_resolution(resolution)
    missing = [p.name for p in tif_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"[{resolution}] Missing TIFFs after extraction: {missing}")

    return tif_paths

# ---------------------------------------------------------------------------
# SAMPLING
# ---------------------------------------------------------------------------

def clean_value(v, nodata=None):
    if v is None:
        return float("nan")
    try:
        v = float(v)
    except Exception:
        return float("nan")
    if nodata is not None and v == nodata:
        return float("nan")
    if v <= -3.0e38:
        return float("nan")
    return v

def sample_raster(tif_path: Path, coordinates: list[tuple[float, float]], method: str) -> list[float]:
    resampling_map = {"nearest": Resampling.nearest, "bilinear": Resampling.bilinear}

    with rasterio.open(tif_path) as src:
        nodata = src.nodata

        if method == "nearest":
            sampled = [x[0] for x in src.sample(coordinates)]
        else:
            with WarpedVRT(src, resampling=resampling_map[method]) as vrt:
                sampled = [x[0] for x in vrt.sample(coordinates)]

    return [clean_value(v, nodata=nodata) for v in sampled]

# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------

def build_output_paths(resolution: str) -> tuple[Path, Path]:
    output_csv = OUTPUT_DIR / f"Anastatica_all_worldclim_bioclim_{resolution}.csv"
    metadata_txt = OUTPUT_DIR / f"Anastatica_all_worldclim_metadata_{resolution}.txt"
    return output_csv, metadata_txt

def write_metadata(resolution: str, source_url: str, output_csv: Path, metadata_txt: Path, n_points: int) -> None:
    text = f"""WorldClim extraction metadata
============================

Input CSV:
{INPUT_CSV}

Output CSV:
{output_csv}

Metadata file:
{metadata_txt}

Cache directory:
{CACHE_DIR}

Resolution:
{resolution}

Sampling method:
{RESAMPLING_METHOD}

Latitude column:
{LATITUDE_COLUMN}

Longitude column:
{LONGITUDE_COLUMN}

Rows processed:
{n_points}

Dataset:
WorldClim v2.1 historical baseline (1970-2000)

Source URL:
{source_url}

Notes:
- Existing ZIP files are reused.
- Existing extracted TIFFs are reused when present.
- No bio_name or bio_detail columns are written.
- BIO1-BIO19 are written as numeric columns only.
"""
    metadata_txt.parent.mkdir(parents=True, exist_ok=True)
    metadata_txt.write_text(text, encoding="utf-8")

def extract_for_resolution(base_df: pd.DataFrame, coordinates: list[tuple[float, float]], resolution: str) -> None:
    log(f"[{resolution}] Starting extraction.")
    tif_paths = ensure_worldclim_data(resolution)

    df = base_df.copy()

    for bio_num, tif_path in zip(BIO_CODES, tif_paths):
        log(f"[{resolution}] Sampling BIO{bio_num} from {tif_path.name}")
        df[f"bio{bio_num}"] = sample_raster(tif_path, coordinates, RESAMPLING_METHOD)

    output_csv, metadata_txt = build_output_paths(resolution)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False, quoting=csv.QUOTE_MINIMAL)

    write_metadata(
        resolution=resolution,
        source_url=WORLDCLIM_URLS[resolution],
        output_csv=output_csv,
        metadata_txt=metadata_txt,
        n_points=len(df),
    )

    log(f"[{resolution}] Wrote: {output_csv}")
    log(f"[{resolution}] Wrote: {metadata_txt}")

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    validate_configuration()

    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    log(f"Reading input: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    df = validate_input(df)

    coordinates = list(zip(df[LONGITUDE_COLUMN].tolist(), df[LATITUDE_COLUMN].tolist()))

    for resolution in WORLDCLIM_URLS:
        extract_for_resolution(df, coordinates, resolution)

    log("All resolutions completed successfully.")

if __name__ == "__main__":
    main()