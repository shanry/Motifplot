import os
import subprocess

import argparse
import sys

import pandas as pd

import RNA
from pdfCropMargins import crop


OUT_DIR = './output_motifs'

COLOR_GLOBAL = "GREEN"


def replace_substring(string, start, end, replacement):
    return string[:start] + replacement + string[end:]


def shrink(i, j, bpairs, SIZE = 2):
    inew, jnew = i, j
    for pair in bpairs:
        delta = pair[1] - pair[0] - 1 - SIZE
        if pair[1] < i:
            inew -= delta
            jnew -= delta
        elif pair[0] > i and pair[1] < j:
            jnew -= delta
        elif pair[0] == i or pair[1] == j:
            jnew -= delta
    return inew, jnew


def extract_pairs_list(ss):
    pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c=='.':
            pass
        elif c=="(":
            stack.append(i)
        elif c==")":
            j = stack.pop()
            pairs.append((j, i))
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs


def find_all_positions(text, substring):
  """Finds all starting indices of a substring in a text.

  Args:
    text: The main text string.
    substring: The substring to search for.

  Returns:
    A list of indices where the substring starts.
  """

  start = 0
  while True:
    start = text.find(substring, start)
    if start == -1:
      break
    yield start
    start += len(substring)


def dotbracket_to_dict(dotbracket):
    """Convert dot-bracket notation to a dictionary with positions as keys and base pairs as values."""
    mstr = dotbracket.replace('(*)', '(***)')
    has_external = False
    if mstr[0] == '5' and mstr[-1] == '3':
        has_external = True
        mstr = mstr[1:-1]
    positions = list(find_all_positions(mstr, '(***)'))
    if has_external:
        bpairs = [(-1, len(mstr))]
    else:
        bpairs = [(0, len(mstr)-1)]
    ipairs = []
    for i in positions:
        bpairs.append((i, i+4))
    mstr = mstr.replace('(***)', '(...)')
    pairs_all = extract_pairs_list(mstr)
    for i, j in pairs_all:
        if (i, j) not in bpairs:
            ipairs.append((i, j))
    motif_dict = dict()
    motif_dict['dot-bracket'] = [dotbracket]
    motif_dict['y_star'] = mstr
    motif_dict['y_sub'] = mstr
    motif_dict['id'] = "short_enum"
    motif_dict['bpairs'] = bpairs
    motif_dict['ipairs'] = ipairs
    motif_dict['length'] = len(mstr)
    motif_dict['cardinality'] = len(ipairs) + 1
    motif_dict['has_external'] = has_external

    return motif_dict


# extend the outmost pair with two unpaired bases
def motif_to_plotstr(m, name=None):
    if name is None:
        m['motif_id'] = "motif"
    else:
        m['motif_id'] = name
    # shrink y
    # ynew = m['y_sub']
    m['bpairs'] = sorted(m['bpairs'])
    m['ipairs'] = sorted(m['ipairs'])
    if m['bpairs'][0][0] >= 0:
        ysub = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        ynew = m['y_star'][m['bpairs'][0][0]: m['bpairs'][0][1]+1]
        start = m['bpairs'][0][0]
    else:
        ysub = m['y_star']
        ynew = m['y_star']
        # m['bpairs'][0] = [0, len(ysub)-1]
        start = 0
    # print(m['bpairs'])
    # print(m['ipairs'])
    # print(ynew)
    SIZE = 3
    delta = 0
    # start = m['bpairs'][0][0]
    for pair in m['bpairs'][1:]:
        # print(pair, pair[1] - pair[0] - 1 - 2)
        # print(len(ynew))
        ynew = replace_substring(ynew, pair[0]+1-start-delta, pair[1]-start-delta, "...")
        # print(len(ynew))
        delta += pair[1] - pair[0] - 1 - SIZE
        # print(f"delta: {delta}")
        # print(f"ydiff: {len(m['y_sub']) - len(ynew)}")
        assert len(ysub) - len(ynew) == delta, f"{len(ysub), len(ynew)}"
        # print()
    # shrink pairs
    bpairs = []
    ipairs = []
    for i, j in m['bpairs']:
        if i >= 0:
            bpairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
        else:
            bpairs.append((i, j))
    for i, j in m['ipairs']:
        ipairs.append(shrink(i, j, m['bpairs'][1:], SIZE))
    i0, j0 = bpairs[0]
    i_base = i0
    flag = False
    if i0 >= 0:
        ynew = "." + ynew + "." # extend the outmost pair with two unpaired bases
        i_base -= 1
        greenplot = f"{i0+1-i_base} {j0+1-i_base} {COLOR_GLOBAL} Fomark "
        flag = True
    else:
        assert len(m['y_star']) == len(m['y_sub'])
        greenplot = f"{1} {len(ynew)} {COLOR_GLOBAL} Fomark "
        i_base = 0
    whiteplot = ""
    for i, j in bpairs[1:]:
        whiteplot += f"{i+1-i_base} {j+1-i_base} WHITE Fomark "
    pairplot = ""
    if bpairs[0][0] >= 0:
        pairplot += f"{i0+1-i_base} {j0+1-i_base} 0.667 0.5 colorpair "
    for i, j in bpairs[1:]:
        pairplot += f"{i+1-i_base} {j+1-i_base} 0.667 0.5 colorpair "
    # for i, j in ipairs:  # internal pairs
    #     pairplot += f"{i+1-i_base} {j+1-i_base} 0.1667 1.0 colorpair "
    poststr = ""
    for i, j in bpairs[1:]:
        poststr += f"{i+1-i_base+1} {j+1-i_base-1} 10 WHITE omark "
    if flag:
        poststr += f"1 1 10 WHITE omark "
        poststr += f"{len(ynew)} {len(ynew)} 10 WHITE omark "
    plotstr = ",".join([m['motif_id'], ynew, greenplot + whiteplot + pairplot, poststr])
    # plotstr = ",".join([m['motif_id'], ynew, greenplot + whiteplot, poststr])  # remove pairplot
    # plotstr = ",".join([m['motif_id'], ynew, whiteplot, poststr])
    plotstr = '"' + plotstr + '"'
    return plotstr


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
                check=True
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
                print(f"Plot saved as {saved_path}")
            else:
                print("No valid output file found in the script output.")
            return True, saved_path, stdout, stderr
        else:
            proc = subprocess.run(
                ["./scripts/draw_motif.sh", plot_str],
                check=True
            )
            return proc.returncode == 0, None, None, None
    except subprocess.CalledProcessError as e:
        return False, None, getattr(e, 'stdout', None), getattr(e, 'stderr', None)


def example():
    db = "((*).(*)(.(*).))"
    motif_dict = dotbracket_to_dict(db)
    print(f"Motif Dictionary: {motif_dict}")

    plotstr = motif_to_plotstr(motif_dict, name="example")
    print(f"Plot String:")
    print(f"{plotstr}")

    draw_motif(plotstr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RNA Motif Plotting Tool')
    parser.add_argument('-e', "--example", action='store_true', help='Run an example motif plot')
    parser.add_argument('-n', '--name', type=str, default='motif', help='Name for the motif')
    parser.add_argument('-o', '--out_dir', type=str, default='./output_motifs', help='Output directory for plots')
    parser.add_argument('--online', action='store_true', help='Read dot-bracket notations from standard input')
    args = parser.parse_args()

    if args.online:
        count = 0
        for line in sys.stdin:
            motif = line.strip()
            if motif:
                motif_dict = dotbracket_to_dict(motif)
                print(f"Motif Dictionary: {motif_dict}")
                plotstr = motif_to_plotstr(motif_dict, name=f"{args.name}_{count}")
                print(f"Plot String:")
                print(f"{plotstr}")
                draw_motif(plotstr)
                count += 1
    elif args.example:
        example()