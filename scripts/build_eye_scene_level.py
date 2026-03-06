#!/usr/bin/env python3
"""Build eye-tracking scene-level summary tables.

Input assumptions (v1):
- One CSV = one subject x one scene.
- Filename contains subject_id and scene_id (customizable via regex).
- CSV contains raw sample/event rows exported by the eye tracker.

This script intentionally does NOT depend on EEG internals. It produces a
scene-level table that can later be merged into EEG batch summaries by
(subject_id, scene_id).

Example:
  python scripts/build_eye_scene_level.py \
    --input /path/to/eye_csvs \
    --output eye_outputs/summary/all_subjects_eye_scene_level.csv

Optional regex example:
  --filename-regex '(?P<subject_id>S\d+).*?(?P<scene_id>\d+)'
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

DEFAULT_FILENAME_REGEX = r"(?P<subject_id>[A-Za-z0-9_-]+).*?(?:scene|round|trial)?[_-]?(?P<scene_id>\d+)"

COLUMN_ALIASES = {
    "tracking_ratio": ["Tracking Ratio[%]"],
    "timestamp_ms": ["Recording Time Stamp[ms]"],
    "event_label": ["Event Label"],
    "annotation": ["Annotation"],
    "trigger_send": ["Triggle Send", "Trigger Send"],
    "trigger_receive": ["Triggle Receive", "Trigger Receive"],
    "validity_left": ["Validity Left"],
    "validity_right": ["Validity Right"],
    "pupil_left_mm": ["Pupil Diameter Left[mm]"],
    "pupil_right_mm": ["Pupil Diameter Right[mm]"],
    "gaze_velocity": ["Gaze Velocity[px/ms]"],
    "fix_idx": ["Fixation Index"],
    "fix_dur_ms": ["Fixation Duration[ms]"],
    "sacc_idx": ["Saccade Index"],
    "sacc_dur_ms": ["Saccade Duration[ms]"],
    "sacc_amp_px": ["Saccade Amplitude[px]"],
    "sacc_vel_avg": ["Saccade Velocity Average[px/ms]"],
    "sacc_vel_peak": ["Saccade Velocity Peak[px/ms]"],
    "blink_idx": ["Blink Index"],
    "blink_dur_ms": ["Blink Duration[ms]"],
}


def find_col(df: pd.DataFrame, aliases: Iterable[str]) -> Optional[str]:
    for a in aliases:
        if a in df.columns:
            return a
    return None


def as_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def unique_count(s: pd.Series) -> float:
    x = pd.to_numeric(s, errors="coerce").dropna()
    return float(x.nunique()) if len(x) else math.nan


def mean_if(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    return float(as_numeric(df[col]).mean())


def max_if(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    return float(as_numeric(df[col]).max())


def unique_if(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    return unique_count(df[col])


def validity_ratio(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    x = as_numeric(df[col])
    if x.dropna().empty:
        return math.nan
    # heuristic: 0 is valid for many exports; if range differs, preserve raw mean
    vals = x.dropna()
    if set(vals.unique()).issubset({0, 1}):
        return float((vals == 0).mean())
    return float(vals.mean())


def parse_filename(path: Path, pattern: re.Pattern[str]) -> Tuple[Optional[str], Optional[int]]:
    m = pattern.search(path.stem)
    if not m:
        return None, None
    sid = m.groupdict().get("subject_id")
    scene = m.groupdict().get("scene_id")
    try:
        scene_id = int(scene) if scene is not None else None
    except Exception:
        scene_id = None
    return sid, scene_id


def summarize_file(path: Path, pattern: re.Pattern[str]) -> Dict[str, object]:
    df = pd.read_csv(path)
    subject_id, scene_id = parse_filename(path, pattern)

    ts_col = find_col(df, COLUMN_ALIASES["timestamp_ms"])
    trig_send_col = find_col(df, COLUMN_ALIASES["trigger_send"])
    trig_recv_col = find_col(df, COLUMN_ALIASES["trigger_receive"])
    event_col = find_col(df, COLUMN_ALIASES["event_label"])
    ann_col = find_col(df, COLUMN_ALIASES["annotation"])

    start_ms = float(as_numeric(df[ts_col]).min()) if ts_col else math.nan
    end_ms = float(as_numeric(df[ts_col]).max()) if ts_col else math.nan

    pupil_left = mean_if(df, "pupil_left_mm")
    pupil_right = mean_if(df, "pupil_right_mm")
    pupil_mean = pd.Series([pupil_left, pupil_right], dtype="float64").mean(skipna=True)

    out = {
        "subject_id": subject_id,
        "scene_id": scene_id,
        "eye_source_file": str(path),
        "eye_n_rows": int(len(df)),
        "eye_event_start_ms": start_ms,
        "eye_event_end_ms": end_ms,
        "eye_duration_ms": (end_ms - start_ms) if pd.notna(start_ms) and pd.notna(end_ms) else math.nan,
        "eye_tracking_ratio": mean_if(df, "tracking_ratio"),
        "eye_valid_left_ratio": validity_ratio(df, "validity_left"),
        "eye_valid_right_ratio": validity_ratio(df, "validity_right"),
        "eye_mean_pupil_left_mm": pupil_left,
        "eye_mean_pupil_right_mm": pupil_right,
        "eye_mean_pupil_mm": float(pupil_mean) if pd.notna(pupil_mean) else math.nan,
        "eye_mean_gaze_velocity_px_ms": mean_if(df, "gaze_velocity"),
        "eye_fix_count": unique_if(df, "fix_idx"),
        "eye_mean_fix_dur_ms": mean_if(df, "fix_dur_ms"),
        "eye_sacc_count": unique_if(df, "sacc_idx"),
        "eye_mean_sacc_dur_ms": mean_if(df, "sacc_dur_ms"),
        "eye_mean_sacc_amp_px": mean_if(df, "sacc_amp_px"),
        "eye_mean_sacc_vel_px_ms": mean_if(df, "sacc_vel_avg"),
        "eye_peak_sacc_vel_px_ms": max_if(df, "sacc_vel_peak"),
        "eye_blink_count": unique_if(df, "blink_idx"),
        "eye_mean_blink_dur_ms": mean_if(df, "blink_dur_ms"),
        "eye_has_event_label": int(event_col is not None),
        "eye_has_annotation": int(ann_col is not None),
        "eye_has_trigger_send": int(trig_send_col is not None),
        "eye_has_trigger_receive": int(trig_recv_col is not None),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="folder containing eye CSV files")
    ap.add_argument("--output", required=True, help="output CSV path")
    ap.add_argument("--filename-regex", default=DEFAULT_FILENAME_REGEX)
    args = ap.parse_args()

    in_dir = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    pattern = re.compile(args.filename_regex, flags=re.IGNORECASE)
    files = sorted(in_dir.rglob("*.csv"))
    rows: List[Dict[str, object]] = []
    for f in files:
        try:
            rows.append(summarize_file(f, pattern))
        except Exception as e:
            rows.append({
                "subject_id": None,
                "scene_id": None,
                "eye_source_file": str(f),
                "eye_parse_error": str(e),
            })

    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
