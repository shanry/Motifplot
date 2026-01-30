import os
import random
import argparse
import subprocess
# from multiprocessing import Pool

import RNA
from pdfCropMargins import crop


OUT_DIR = './output_structures'


def plot_example():
    y = '(((((((((....)))))))))......((((((................((((((....))))))................))))))..'
    x = 'o' * len(y)
    RNA.PS_rna_plot_a(x, y, 'rna.ps', pre='1 15 8 RED omark', post='20 27 8 GREEN omark')
    # ps2pdf -dEPSCrop rna.ps
    os.system('ps2pdf -dEPSCrop rna.ps')
    # remove the ps file
    os.remove('rna.ps')
    print('Plot saved as rna.pdf')


def plot_structure(y, name=None, out_dir=None):
    if out_dir is None:
        out_dir = OUT_DIR
    x = 'o' * len(y)
    if name is None:
        filename = 'rna' + str(random.randint(0, 9999999))
    else:
        filename = name
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    path_ps = os.path.join(out_dir, filename + '.ps')
    path_pdf = os.path.join(out_dir, filename + '.pdf')
    path_pdf_cropped = os.path.join(out_dir, filename + '_cropped.pdf')
    RNA.PS_rna_plot_a(x, y, path_ps, pre='', post='')
    # Convert PS to PDF
    os.system(f'ps2pdf -dEPSCrop {path_ps} {path_pdf}')
    # Remove the PS file
    os.remove(path_ps)

    crop(["-u", "-s", f"{path_pdf}", "-o", f"{path_pdf_cropped}"])
    # mv cropped pdf to pdf
    os.rename(path_pdf_cropped, path_pdf)

    print(f'Plot saved as {path_pdf}')


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


def batch_plot(path_jsonl):
    import json
    for i, line in enumerate(open(path_jsonl)):
        data = json.loads(line)
        structure = data['target_structure']
        ID = data['id']
        name = f'{ID}'
        print(f'Plotting {name}: {structure}')
        # plot_structure(structure, name=name)
        plotstr = f'"{name},{structure}"'
        print(f"Plot String:")
        print(f"{plotstr}")
        draw_motif(plotstr, out_dir=OUT_DIR)


def batch_plot_parallel(path_jsonl, num_processes=None):
    import json
    from multiprocessing import Pool
    if num_processes is None:
        num_processes = os.cpu_count()

    tasks = []
    for i, line in enumerate(open(path_jsonl)):
        data = json.loads(line)
        structure = data['target_structure']
        ID = data['id']
        name = f'{ID}'
        # tasks.append((structure, name))
        plotstr = f'"{name},{structure}"'
        tasks.append((plotstr, True, OUT_DIR))

    with Pool(processes=num_processes) as pool:
        pool.starmap(draw_motif, tasks)

    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RNA Structure Plotting Tool')
    # path
    parser.add_argument('-p', '--path', type=str, default='structures.jsonl', help='Path to the JSONL file for batch plotting')
    parser.add_argument('-l', '--layout', type=int, default=0, help='Layout type (0: SIMPLE, 1: NAVIEW, 2: CIRCULAR, 3: TURTLE, 4: PUZZLER)')
    parser.add_argument('-e', '--example', action='store_true', help='Run example RNA structure plot')
    parser.add_argument('-y', '--structure', type=str, help='RNA structure in dot-bracket notation')
    parser.add_argument('-n', '--name', type=str, default="y", help='Name for the output file')
    parser.add_argument('-r', '--parallel', action='store_true', help='Run batch plotting in parallel')
    parser.add_argument('-o', '--out_dir', type=str, default='./output_structures', help='Output directory for plots')


    args = parser.parse_args()

    RNA.cvar.rna_plot_type = args.layout
    print(f'Using layout type: {RNA.cvar.rna_plot_type}')

    OUT_DIR = args.out_dir
    print(f'Output directory: {OUT_DIR}')

    if args.example:
        plot_example()
    elif args.structure:
        # plot_structure(args.structure, args.name)
        plotstr = f'"{args.name},{args.structure}"'
        print(f"Plot String:")
        print(f"{plotstr}")
        success, saved_path, stdout, stderr = draw_motif(plotstr, out_dir=OUT_DIR)
    else:    
        if args.path:
            if args.parallel:
                print(f'Running batch plotting in parallel using {os.cpu_count()} processes.')
                batch_plot_parallel(args.path)
            else:
                print('Running batch plotting in sequential mode.')
                batch_plot(args.path)
        else:
            print("No structure provided. Use -y to specify a structure or -p for batch plotting.")