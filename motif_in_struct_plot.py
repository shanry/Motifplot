import os
import subprocess
import argparse

import pandas as pd

import RNA
from pdfCropMargins import crop

import json

OUT_DIR = "./output_mins"  # miny: m (moitf) in s (structure)

# a list of colors for plotting motifs
COLORS = [ "GREEN", "1 1 0", "0 1 1", "1 0.65 0", "1 0 1"]


# plot motifs in original structure
def get_prestr(
    y_star, bpairs, ipairs=None, color="GREEN"
):  # bpairs: boundary pairs, ipairs: internal pairs
    i0, j0 = bpairs[0]
    if i0 >= 0:
        greenplot = f"{i0+1} {j0+1} {color} Fomark "
    else:
        greenplot = f"{1} {len(y_star)} {color} Fomark "
        i0 = 0
        j0 = len(y_star) - 1
    whiteplot = ""
    for i, j in bpairs[1:]:
        whiteplot += f"{i+1} {j+1} WHITE Fomark "
    pairplot = ""
    if bpairs[0][0] >= 0:
        pairplot += f"{i0+1} {j0+1} 0.667 0.5 colorpair "
    for i, j in bpairs[1:]:
        pairplot += f"{i+1} {j+1} 0.667 0.5 colorpair "
    if ipairs is not None:
        for i, j in ipairs:
            pairplot += f"{i+1} {j+1} 0.1667 1.0 colorpair "
    # plotstr = ",".join([m['ymotif_id'], m['y_star'], greenplot + whiteplot + pairplot])
    prestr = greenplot + whiteplot + pairplot
    return prestr


def draw_motif(plot_str, capture=True, out_dir=None):
    plot_str = plot_str.strip('"')
    if out_dir is None:
        out_dir = OUT_DIR
    try:
        if capture:
            proc = subprocess.run(
                ["./scripts/draw_motif.sh", plot_str],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
            stdout = proc.stdout
            stderr = proc.stderr
            last_line = stdout.rstrip().splitlines()[-1] if stdout.strip() else ""
            saved_path = None
            if last_line and os.path.exists(last_line):
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir, exist_ok=True)
                basename = os.path.basename(last_line)
                output_path = os.path.join(out_dir, basename)
                os.replace(last_line, output_path)
                saved_path = output_path
                # print(f"Plot saved as {saved_path}")
            else:
                print("No valid output file found in the script output.")
            return True, saved_path, stdout, stderr
        else:
            proc = subprocess.run(["./scripts/draw_motif.sh", plot_str], check=True)
            return proc.returncode == 0, None, None, None
    except subprocess.CalledProcessError as e:
        return False, None, getattr(e, "stdout", None), getattr(e, "stderr", None)


def example():
    y_star = ".....((((((((((((....))))((((....))))))))((((((((....))))((((....))))))))))))...................."
    bpairs_list = [[(8, 73), (11, 38), (41, 72)], [(12, 37), (15, 22), (25, 36)], [(44, 69), (47, 54), (57, 68)]]
    for i, bpairs in enumerate(bpairs_list):
        prestr = get_prestr(y_star, bpairs)
        ID = f"example_{i}"
        plotstr = f'"{ID},{y_star},{prestr}"'
        if args.verbose:
            print(f"Plot String:")
            print(f"{plotstr}")
        success, saved_path, _, _ = draw_motif(
            plotstr, capture=True, out_dir=OUT_DIR
        )
        if success:
            print(f"Plot saved as {saved_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA Motif Plotting Tool")
    parser.add_argument('-p', '--path', type=str, default='', help='Path to the file containing dot-bracket notations')
    parser.add_argument('-d', '--dotbracket', type=str, default='', help='RNA structure in dot-bracket notation')
    parser.add_argument('-e', "--example", action='store_true', help='Run an example motif plot')
    parser.add_argument('-n', '--name', type=str, default='motif', help='Name for the motif')
    parser.add_argument('--csv', action='store_true', help='Indicates that the input file is a CSV file with a "structure" column')
    parser.add_argument('-o', '--out_dir', type=str, default='./output_mins', help='Output directory for plots')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()

    OUT_DIR = args.out_dir
    print(f"Output directory: {OUT_DIR}")

    example()