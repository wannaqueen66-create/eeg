#!/usr/bin/env python3
"""
Rebuild Task8 (group regulation modes) numeric tables + PNG cards from existing QC outputs.

Usage:
  python3 scripts/rebuild_task8_numeric.py \
    --root /root/analysis-2-EEG \
    --out  /root/analysis-2-EEG/task8_group_regulation_modes
"""

import argparse
import csv
import os
import sys

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception as e:
    print("[ERROR] Pillow is required for PNG export. Install with: python3 -m pip install --break-system-packages pillow")
    raise


def read_csv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def find(rows, key, val):
    for r in rows:
        if r.get(key) == val:
            return r
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Analysis root containing task3/task4 QC tables")
    ap.add_argument("--out", required=True, help="Output folder for task8_group_regulation_modes")
    args = ap.parse_args()

    root = args.root
    out = args.out
    os.makedirs(f"{out}/tables/qc", exist_ok=True)
    os.makedirs(f"{out}/figures/qc", exist_ok=True)
    os.makedirs(f"{out}/reports/qc", exist_ok=True)

    analyses = ["experience", "sportfreq"]
    metrics = ["O_theta", "F_theta", "O_alpha", "O_beta"]

    rows_out = []
    for a in analyses:
        for m in metrics:
            r3 = read_csv(f"{root}/task3_trialindex_lmm/tables/qc/{a}/lmm_fixed_effects_{m}_qc.csv")
            b_trial = float(find(r3, "Name", "TrialIndex")["Estimate"])
            r_int = find(r3, "Name", "Group_High:TrialIndex")
            b_gxt = float(r_int["Estimate"])
            p_gxt = float(r_int["pValue"])
            low_slope = b_trial
            high_slope = b_trial + b_gxt
            slope_gap = high_slope - low_slope

            r4 = read_csv(f"{root}/task4_core_lmm_suite/factor_WWR/tables/qc/{a}/model2_two_way_{m}_qc_fixed_effects.csv")
            b_c = float(find(r4, "Name", "Complexity_ComplexityHigh")["Estimate"])
            r_cg = find(r4, "Name", "Complexity_ComplexityHigh:Group_High")
            b_cg = float(r_cg["Estimate"])
            p_cg = float(r_cg["pValue"])
            low_comp = b_c
            high_comp = b_c + b_cg
            comp_gap = high_comp - low_comp

            r41 = read_csv(f"{root}/task4_core_lmm_suite/factor_WWR/tables/qc/{a}/model1_main_effects_{m}_qc_fixed_effects.csv")
            rg = find(r41, "Name", "Group_High")
            b_group = float(rg["Estimate"])
            p_group = float(rg["pValue"])

            rows_out.append(
                {
                    "analysis": a,
                    "metric": m,
                    "low_trial_slope": low_slope,
                    "high_trial_slope": high_slope,
                    "slope_gap_high_minus_low": slope_gap,
                    "p_group_by_trial": p_gxt,
                    "low_complexity_effect": low_comp,
                    "high_complexity_effect": high_comp,
                    "complexity_effect_gap_high_minus_low": comp_gap,
                    "p_complexity_by_group": p_cg,
                    "group_main_diff_controlled_high_minus_low": b_group,
                    "p_group_main_controlled": p_group,
                }
            )

    csv_path = f"{out}/tables/qc/group_regulation_numeric_effects_qc.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)

    # Summary PNG
    W, H = 2600, 900
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    font = ImageFont.load_default()

    d.text((20, 15), "Task8 QC - Group Regulation Numeric Effects (High vs Low)", fill="black", font=font)
    d.text((20, 35), "Columns: slope gap (Task3), complexity modulation gap (Task4 model2), group main diff controlled (Task4 model1)", fill="black", font=font)

    headers = [
        "analysis",
        "metric",
        "slope_low",
        "slope_high",
        "gap_slope(H-L)",
        "p(GxTrial)",
        "compEff_low",
        "compEff_high",
        "gap_comp(H-L)",
        "p(CxG)",
        "group_main(H-L)",
        "p(Group)",
    ]
    colx = [20, 150, 280, 410, 560, 730, 900, 1060, 1230, 1410, 1560, 1770]
    for x, h in zip(colx, headers):
        d.text((x, 70), h, fill="navy", font=font)

    y = 95
    for r in rows_out:
        vals = [
            r["analysis"],
            r["metric"],
            f"{r['low_trial_slope']:.6f}",
            f"{r['high_trial_slope']:.6f}",
            f"{r['slope_gap_high_minus_low']:.6f}",
            f"{r['p_group_by_trial']:.4g}",
            f"{r['low_complexity_effect']:.6f}",
            f"{r['high_complexity_effect']:.6f}",
            f"{r['complexity_effect_gap_high_minus_low']:.6f}",
            f"{r['p_complexity_by_group']:.4g}",
            f"{r['group_main_diff_controlled_high_minus_low']:.6f}",
            f"{r['p_group_main_controlled']:.4g}",
        ]
        if r["p_group_by_trial"] < 0.05 or r["p_complexity_by_group"] < 0.05:
            d.rectangle((15, y - 2, 1900, y + 12), outline="red", width=1)
        for x, v in zip(colx, vals):
            d.text((x, y), v, fill="black", font=font)
        y += 24

    scale = 18000
    by = 100
    d.text((1950, 70), "Visual bars (signed effect sizes)", fill="navy", font=font)
    for r in rows_out:
        d.text((1950, by), f"{r['analysis']} {r['metric']}", fill="black", font=font)
        x0 = 2200
        d.line((x0 - 120, by + 8, x0 + 120, by + 8), fill="gray", width=1)
        for i, (val, color) in enumerate(
            [
                (r["slope_gap_high_minus_low"], "blue"),
                (r["complexity_effect_gap_high_minus_low"], "green"),
                (r["group_main_diff_controlled_high_minus_low"], "purple"),
            ]
        ):
            y0 = by + 2 + i * 6
            y1 = y0 + 4
            dx = int(val * scale)
            if dx >= 0:
                d.rectangle((x0, y0, x0 + dx, y1), fill=color)
            else:
                d.rectangle((x0 + dx, y0, x0, y1), fill=color)
        by += 85

    d.text((1950, 780), "blue=slope gap(H-L); green=complexity-effect gap(H-L); purple=group main diff(H-L)", fill="black", font=font)

    png_path = f"{out}/figures/qc/task8_group_regulation_numeric_summary_qc.png"
    img.save(png_path)

    for r in rows_out:
        card = Image.new("RGB", (1200, 520), "white")
        dr = ImageDraw.Draw(card)
        dr.text((20, 20), f"Task8 QC Numeric Card | {r['analysis']} | {r['metric']}", fill="black", font=font)
        dr.text((20, 55), f"Trial slope (Low): {r['low_trial_slope']:.6f}", fill="black", font=font)
        dr.text((20, 75), f"Trial slope (High): {r['high_trial_slope']:.6f}", fill="black", font=font)
        dr.text((20, 95), f"Slope gap (High-Low): {r['slope_gap_high_minus_low']:.6f} | p(GxTrial)={r['p_group_by_trial']:.4g}", fill="black", font=font)
        dr.text((20, 135), f"Complexity effect in Low: {r['low_complexity_effect']:.6f}", fill="black", font=font)
        dr.text((20, 155), f"Complexity effect in High: {r['high_complexity_effect']:.6f}", fill="black", font=font)
        dr.text((20, 175), f"Complexity-effect gap (High-Low): {r['complexity_effect_gap_high_minus_low']:.6f} | p(CxG)={r['p_complexity_by_group']:.4g}", fill="black", font=font)
        dr.text((20, 215), f"Group main diff controlled (High-Low): {r['group_main_diff_controlled_high_minus_low']:.6f} | p(Group)={r['p_group_main_controlled']:.4g}", fill="black", font=font)

        cx = 760
        scale2 = 24000
        dr.text((620, 55), "Signed effect bars", fill="navy", font=font)
        dr.line((cx - 250, 110, cx + 250, 110), fill="gray", width=1)

        def bar(yc, val, color, label):
            dx = int(val * scale2)
            if dx >= 0:
                dr.rectangle((cx, yc - 8, cx + dx, yc + 8), fill=color)
            else:
                dr.rectangle((cx + dx, yc - 8, cx, yc + 8), fill=color)
            dr.text((20, yc - 8), label, fill="black", font=font)

        bar(130, r["slope_gap_high_minus_low"], "blue", "Slope gap (H-L)")
        bar(180, r["complexity_effect_gap_high_minus_low"], "green", "Complexity gap (H-L)")
        bar(230, r["group_main_diff_controlled_high_minus_low"], "purple", "Group main diff (H-L)")

        flag = "sig" if (r["p_group_by_trial"] < 0.05 or r["p_complexity_by_group"] < 0.05) else "ns"
        dr.text((20, 460), f"Interaction evidence flag: {flag}", fill="red" if flag == "sig" else "black", font=font)

        outp = f"{out}/figures/qc/task8_numeric_card_qc_{r['analysis']}_{r['metric']}.png"
        card.save(outp)

    # tiny report
    with open(f"{out}/reports/qc/group_regulation_numeric_report_qc.md", "w", encoding="utf-8") as f:
        f.write("# Task8 Numeric Rebuild (QC)\n\n")
        f.write(f"- Source root: `{root}`\n")
        f.write(f"- Output root: `{out}`\n")
        f.write("- Wrote CSV: `tables/qc/group_regulation_numeric_effects_qc.csv`\n")
        f.write("- Wrote summary PNG: `figures/qc/task8_group_regulation_numeric_summary_qc.png`\n")
        f.write("- Wrote 8 cards: `figures/qc/task8_numeric_card_qc_<analysis>_<metric>.png`\n")

    print("[OK] Rebuilt Task8 numeric outputs")


if __name__ == "__main__":
    main()
