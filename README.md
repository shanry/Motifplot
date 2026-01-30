# MotifPlot
RNA Visualization

## Installation
1.  **Clone the repository:**
```bash
git clone https://github.com/shanry/Motifplot.git
```
2.  **Install dependencies:**: This project is developed using Python 3.11, and should also work for other Python versions. The ViennaRNA Package (2.7.0) is needed: https://www.tbi.univie.ac.at/RNA/

```bash
pip install -r requirements.txt
```
```bash
export VIENNA=path/to/vrna270/
```

## Plot Structure
```
python struct_plot.py --structure "(((((((....(((...........)))((((((((..(((((((((((((((((((...(((((......))))).)))))).)))))))))))))..))))))))..)))))))"
```

## Plot Motif
```
echo "(.(.((.(*)))).)" | python motif_plot.py --online
```

## Plot Motif in Structure
```
python motif_in_struct_plot.py
```