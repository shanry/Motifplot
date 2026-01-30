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
    y_star = "........(((((((((((....))))..((((....))))..((((....)))))))))))..(((((((...))))..((((...))))..((((...)))))))........."
    motif_bpairs_list = [
        (
            "(((*))..(*)..(*))",
            [
                [(14, 55), (16, 25), (29, 40), (43, 54)],
                [(66, 104), (68, 76), (80, 90), (93, 103)],
            ],
        ),
        (
            "((*)..((*))..(*))",
            [
                [(14, 55), (15, 26), (30, 39), (43, 54)],
                [(66, 104), (67, 77), (81, 89), (93, 103)],
            ],
        ),
        (
            "((*)..(*)..((*)))",
            [
                [(14, 55), (15, 26), (29, 40), (44, 53)],
                [(66, 104), (67, 77), (80, 90), (94, 102)],
            ],
        ),
        (
            "(((*)..(*)..(*)))",
            [
                [(13, 56), (15, 26), (29, 40), (43, 54)],
                [(65, 105), (67, 77), (80, 90), (93, 103)],
            ],
        )
    ]
    y_star = ".....((((((((((((....))))((((....))))))))((((((((....))))((((....))))))))))))...................."
    motif_bpairs_list = [
        (
            "((((*)))(*))",
            [[(8, 73), (11, 38), (41, 72)], [(12, 37), (15, 22), (25, 36)], [(44, 69), (47, 54), (57, 68)]]
        )
    ]
    n_ID = 0
    for _, bpairs_list in motif_bpairs_list:
        print(f"Plotting motif with bpairs: {bpairs_list}")
        for bpairs in bpairs_list:
            prestr = get_prestr(y_star, bpairs)
            ID = f"example_{n_ID}"
            n_ID += 1
            plotstr = f'"{ID},{y_star},{prestr}"'
            if args.verbose:
                print(f"Plot String:")
                print(f"{plotstr}")
            success, saved_path, _, _ = draw_motif(
                plotstr, capture=True, out_dir=OUT_DIR
            )
            if success:
                print(f"Plot saved as {saved_path}")


def example_multi():
    eterna_57 = {"motifs":[[[-1,36],[3,25]],[[3,25],[6,22]],[[6,22]]],"structure":"((.(..(.(....).(....).)..).(....).))"}
    eterna_57_2 = {"motifs":[[[-1,36],[1,34]],[[1,34],[6,22],[27,32]],[[6,22]],[[27,32]]],"structure":"((.(..(.(....).(....).)..).(....).))"}
    eterna_57_3 = {"motifs":[[[-1,36],[3,25],[27,32]],[[3,25],[8,13],[15,20]],[[27,32]],[[8,13]],[[15,20]]],"structure":"((.(..(.(....).(....).)..).(....).))"}
    twomotifs = {"motifs":[[[-1,20],[5,12]],[[5,12]]],"structure":"(((..((....))...).))"}
    eterna_52 = {"motifs":[[[-1,80],[9,18],[33,42],[49,58]],[[9,18]],[[33,42]],[[49,58]]],"structure":".((......((......))......((......((......))......((......))......))......))....."}
    eterna_87 = {"motifs":[[[-1,97],[6,90]],[[6,90],[8,88]],[[8,88],[16,84]],[[16,84],[48,55],[59,66]],[[48,55]],[[59,66]]],"structure":"....((((((((....((..(...................)..((...((....))...((....))........))......))))))))))...."}
    eterna_chickenfeet = {"motifs":[[[-1,67],[8,30],[37,59]],[[8,30],[9,29]],[[9,29]],[[37,59],[38,58]],[[38,58]]],"structure":".....((.((.((...)).((...))...)))).((.((.((...)).((...))...))))....."}
    eterna_50 = {"motifs":[[[-1,105],[5,44],[54,66],[68,77]],[[5,44],[6,43]],[[6,43],[11,26],[30,37]],[[11,26],[12,25]],[[12,25]],[[30,37]],[[54,66]],[[68,77]]],"structure":".....((((.(((((....))....))).(((....))).)).))....((((.(((....)...)).(((....))).))..))...................."}
    eterna_60 = {"motifs":[[[-1,105],[7,82]],[[7,82],[16,73]],[[16,73],[18,71]],[[18,71],[24,62]],[[24,62],[31,56]],[[31,56],[33,54]],[[33,54],[39,48]],[[39,48]]],"structure":".....((((.(...(.((((.(.((((.((.((((.((.(((....))).).)))))).)))))....).)))).)...).))))...................."}
    eterna_67 = {"motifs":[[[-1,61],[19,25]],[[19,25]]],"structure":"......(.........((((.....)))).........)......................"}
    eterna_72 = {"motifs":[[[-1,73],[6,66]],[[6,66],[8,64]],[[8,64],[19,34],[36,53]],[[19,34],[20,33]],[[20,33]],[[36,53],[38,51]],[[38,51]]],"structure":"..((((((((.......(.((((((....)))))).(((((((....))))))).).......)))))))).."}
    decomposition_list = [eterna_57,]# eterna_57_2, eterna_57_3, twomotifs, eterna_52, eterna_87, eterna_chickenfeet, eterna_50,
                          #eterna_60, eterna_67]
    n_ID = 0
    for i_ID, decomposition in enumerate(decomposition_list):
        if i_ID < len(decomposition_list) - 1:
            continue
        y_star = decomposition["structure"]
        bpairs_list = decomposition["motifs"]
        prestr_all = ""
        for ib, bpairs in enumerate(bpairs_list):
            color = COLORS[ib % len(COLORS)]
            prestr = get_prestr(y_star, bpairs, color=color)
            prestr_all += f"{prestr}"
        ID = f"decomposition_{i_ID}"
        n_ID += 1
        plotstr = f'"{ID},{y_star},{prestr_all}"'
        if args.verbose:
            print(f"Plot String:")
            print(f"{plotstr}")
        success, saved_path, _, _ = draw_motif(
            plotstr, capture=True, out_dir=OUT_DIR
        )
        if success:
            print(f"Plot saved as {saved_path}")
    return


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
    # example_multi()