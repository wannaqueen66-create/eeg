#!/usr/bin/env python3
"""Build eye-tracking scene-level summary tables.

Project-specific assumptions (current experiment):
- One CSV = one subject x one scene x one *view* segment.
- Scene identity is parsed from the parent folder name, e.g.:
    1-1-1_组1-C1W45
  meaning:
    order_id=1, round_id=1, scene_index_in_round=1
- Subject identity is parsed from the CSV filename and kept as Chinese name
  to match EEG `.set` filenames.
- Marker cutting is NOT needed on the eye-tracking side in v1.

Example:
  python scripts/build_eye_scene_level.py \
    --input /path/to/eye_root \
    --output eye_outputs/summary/all_subjects_eye_scene_level.csv
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

SCENE_DIR_REGEX = re.compile(
    r"(?P<order_id>\d+)-(?P<round_id>\d+)-(?P<scene_index_in_round>\d+)(?:[_-](?P<scene_label>.*))?$"
)

SCENE_LABEL_REGEX = re.compile(r"C(?P<complexity>[01])W(?P<wwr>15|45|75)", flags=re.IGNORECASE)

# raw_高雅清_260201191617_0207215457.csv -> 高雅清
SUBJECT_FROM_FILE_REGEX = re.compile(r"^(?:raw[_-])?(?P<subject_id>[^_]+?)(?:_\d.*)?$")

COLUMN_ALIASES = {
    "record_name": ["Record Name"],
    "user": ["User"],
    "tracking_ratio": ["Tracking Ratio[%]"],
    "timestamp_ms": ["Recording Time Stamp[ms]"],
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
    "sex": ["性别"],
    "lens_right": ["右眼镜片度数"],
    "lens_left": ["左眼镜片度数"],
}


def find_col(df: pd.DataFrame, aliases: Iterable[str]) -> Optional[str]:
    for a in aliases:
        if a in df.columns:
            return a
    return None


def as_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


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


def unique_count_numeric(s: pd.Series) -> float:
    x = pd.to_numeric(s, errors="coerce").dropna()
    return float(x.nunique()) if len(x) else math.nan


def unique_if(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    return unique_count_numeric(df[col])


def first_text_if(df: pd.DataFrame, key: str) -> str:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return ""
    s = df[col].dropna().astype(str)
    if s.empty:
        return ""
    return s.iloc[0].strip()


def validity_ratio(df: pd.DataFrame, key: str) -> float:
    col = find_col(df, COLUMN_ALIASES.get(key, []))
    if not col:
        return math.nan
    x = as_numeric(df[col]).dropna()
    if x.empty:
        return math.nan
    # heuristic: 0/1 often means valid/invalid; otherwise preserve raw mean
    vals = set(x.unique().tolist())
    if vals.issubset({0, 1}):
        return float((x == 0).mean())
    return float(x.mean())


def parse_subject_id(path: Path) -> str:
    stem = path.stem.strip()
    m = SUBJECT_FROM_FILE_REGEX.match(stem)
    if m:
        return m.group("subject_id").strip()
    return stem


def parse_scene_dir(path: Path) -> Dict[str, object]:
    folder = path.parent.name.strip()
    out: Dict[str, object] = {
        "scene_folder": folder,
        "order_id": math.nan,
        "round_id": math.nan,
        "scene_index_in_round": math.nan,
        "scene_label": "",
        "WWR": math.nan,
        "Complexity": math.nan,
    }

    m = SCENE_DIR_REGEX.match(folder)
    if not m:
        return out

    out["order_id"] = int(m.group("order_id"))
    out["round_id"] = int(m.group("round_id"))
    out["scene_index_in_round"] = int(m.group("scene_index_in_round"))
    label = (m.group("scene_label") or "").strip()
    out["scene_label"] = label

    # global scene_id consistent with existing EEG design convention: (round-1)*6 + pos
    out["scene_id"] = int((out["round_id"] - 1) * 6 + out["scene_index_in_round"])

    m2 = SCENE_LABEL_REGEX.search(label)
    if m2:
        out["Complexity"] = int(m2.group("complexity"))
        out["WWR"] = int(m2.group("wwr"))
    return out


def summarize_file(path: Path) -> Dict[str, object]:
    df = pd.read_csv(path)
    subj = parse_subject_id(path)
    scene_meta = parse_scene_dir(path)

    ts_col = find_col(df, COLUMN_ALIASES["timestamp_ms"])
    start_ms = float(as_numeric(df[ts_col]).min()) if ts_col else math.nan
    end_ms = float(as_numeric(df[ts_col]).max()) if ts_col else math.nan

    pupil_left = mean_if(df, "pupil_left_mm")
    pupil_right = mean_if(df, "pupil_right_mm")
    pupil_mean = pd.Series([pupil_left, pupil_right], dtype="float64").mean(skipna=True)

    out = {
        "subject_id": subj,
        "subject_name": subj,
        "eye_user": first_text_if(df, "user"),
        "eye_record_name": first_text_if(df, "record_name"),
        **scene_meta,
        "eye_source_file": str(path),
        "eye_n_rows": int(len(df)),
        "eye_view_start_ms": start_ms,
        "eye_view_end_ms": end_ms,
        "eye_view_duration_ms": (end_ms - start_ms) if pd.notna(start_ms) and pd.notna(end_ms) else math.nan,
        "eye_tracking_ratio": mean_if(df, "tracking_ratio"),
        "eye_valid_left_ratio": validity_ratio(df, "validity_left"),
        "eye_valid_right_ratio": validity_ratio(df, "validity_right"),
        "eye_mean_pupil_left_mm": pupil_left,
        "eye_mean_pupil_right_mm": pupil_right,
        "eye_mean_pupil_mm": float(pupil_mean) if pd.notna(pupil_mean) else math.nan,
        "eye_view_mean_gaze_velocity_px_ms": mean_if(df, "gaze_velocity"),
        "eye_view_fix_count": unique_if(df, "fix_idx"),
        "eye_view_mean_fix_dur_ms": mean_if(df, "fix_dur_ms"),
        "eye_view_sacc_count": unique_if(df, "sacc_idx"),
        "eye_view_mean_sacc_dur_ms": mean_if(df, "sacc_dur_ms"),
        "eye_view_mean_sacc_amp_px": mean_if(df, "sacc_amp_px"),
        "eye_view_mean_sacc_vel_px_ms": mean_if(df, "sacc_vel_avg"),
        "eye_view_peak_sacc_vel_px_ms": max_if(df, "sacc_vel_peak"),
        "eye_view_blink_count": unique_if(df, "blink_idx"),
        "eye_view_mean_blink_dur_ms": mean_if(df, "blink_dur_ms"),
        "eye_sex": first_text_if(df, "sex"),
        "eye_lens_right": first_text_if(df, "lens_right"),
        "eye_lens_left": first_text_if(df, "lens_left"),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="root folder containing scene subfolders of eye CSVs")
    ap.add_argument("--output", required=True, help="output CSV path")
    args = ap.parse_args()

    in_dir = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(in_dir.rglob("*.csv"))
    rows: List[Dict[str, object]] = []
    for f in files:
        try:
            rows.append(summarize_file(f))
        except Exception as e:
            rows.append({
                "subject_id": "",
                "subject_name": "",
                "scene_id": math.nan,
                "scene_folder": f.parent.name,
                "eye_source_file": str(f),
                "eye_parse_error": str(e),
            })

    df_out = pd.DataFrame(rows)
    preferred_cols = [
        "subject_id", "subject_name", "scene_id", "order_id", "round_id", "scene_index_in_round",
        "scene_folder", "scene_label", "WWR", "Complexity",
        "eye_source_file", "eye_n_rows", "eye_view_start_ms", "eye_view_end_ms", "eye_view_duration_ms",
        "eye_tracking_ratio", "eye_valid_left_ratio", "eye_valid_right_ratio",
        "eye_mean_pupil_left_mm", "eye_mean_pupil_right_mm", "eye_mean_pupil_mm",
        "eye_view_mean_gaze_velocity_px_ms", "eye_view_fix_count", "eye_view_mean_fix_dur_ms",
        "eye_view_sacc_count", "eye_view_mean_sacc_dur_ms", "eye_view_mean_sacc_amp_px",
        "eye_view_mean_sacc_vel_px_ms", "eye_view_peak_sacc_vel_px_ms",
        "eye_view_blink_count", "eye_view_mean_blink_dur_ms",
        "eye_user", "eye_record_name", "eye_sex", "eye_lens_right", "eye_lens_left",
        "eye_parse_error",
    ]
    cols = [c for c in preferred_cols if c in df_out.columns] + [c for c in df_out.columns if c not in preferred_cols]
    df_out = df_out.loc[:, cols]
    df_out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
