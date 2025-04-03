import json, os


def ionic_radius(el, state, coordination, radius_type = 'r_ionic'):
    """
    https://github.com/prtkm/ionic-radii
    R. D. Shannon, Revised Effective Ionic Radii and Systematic Studies of Interatomic Distances in Halides and Chalcogenides,
    Acta Crystallographica A32 (1976) 751-767.

    INPUT:
        el (str) - name of element
        state (int) - oxidation state, -1,1,2,3 ...
        coordination (str) - roman numbers, VI, IV, and spin state for transition metals, e.g. VIHS or VILS
        spin (str) 
            'HS' - high spin
            'LS' - low spin 
        radius_type (str) - radius type and two additional fields
            'r_ionic'
            'r_crystal'
            'spin' - check possible spin
            'remark'



    """
    # print(os.path.dirname(__file__))
    with open(os.path.dirname(__file__)+"/shannon-radii.json") as f:
        out = f.read()

    d = json.loads(out)

    # Enter Element, Charge, Coordination and one of - r_crystal, r_ionic, spin, remark

    return d[el][str(state)][coordination][radius_type]

# ionization_energies = { 'D': 13.60, 'OH': 100,
#     'H': 13.60, 'He': 24.59, 'Li': 5.39, 'Be': 9.32, 'B': 8.30, 
#     'C': 11.26, 'N': 14.53, 'O': 13.62, 'F': 17.42, 'Ne': 21.56, 
#     'Na': 5.14, 'Mg': 7.65, 'Al': 5.99, 'Si': 8.15, 'P': 10.49, 
#     'S': 10.36, 'Cl': 12.97, 'Ar': 15.76, 'K': 4.34, 'Ca': 6.11, 
#     'Sc': 6.56, 'Ti': 6.83, 'V': 6.75, 'Cr': 6.77, 'Mn': 7.43, 
#     'Fe': 7.90, 'Co': 7.88, 'Ni': 7.64, 'Cu': 7.73, 'Zn': 9.39, 
#     'Ga': 6.00, 'Ge': 7.90, 'As': 9.79, 'Se': 9.75, 'Br': 11.81, 
#     'Kr': 14.00, 'Rb': 4.18, 'Sr': 5.69, 'Y': 6.22, 'Zr': 6.63, 
#     'Nb': 6.76, 'Mo': 7.09, 'Tc': 7.28, 'Ru': 7.36, 'Rh': 7.46, 
#     'Pd': 8.34, 'Ag': 7.58, 'Cd': 8.99, 'In': 5.79, 'Sn': 7.34, 
#     'Sb': 8.61, 'Te': 9.01, 'I': 10.45, 'Xe': 12.13, 'Cs': 3.89, 
#     'Ba': 5.21, 'La': 5.58, 'Ce': 5.54, 'Pr': 5.47, 'Nd': 5.53, 
#     'Pm': 5.58, 'Sm': 5.64, 'Eu': 5.67, 'Gd': 6.15, 'Tb': 5.86, 
#     'Dy': 5.94, 'Ho': 6.02, 'Er': 6.11, 'Tm': 6.18, 'Yb': 6.25, 
#     'Lu': 5.43, 'Hf': 6.83, 'Ta': 7.55, 'W': 7.86, 'Re': 7.83, 
#     'Os': 8.44, 'Ir': 8.97, 'Pt': 8.96, 'Au': 9.23, 'Hg': 10.44, 
#     'Tl': 6.11, 'Pb': 7.42, 'Bi': 7.29, 'Po': 8.42, 'At': 9.30, 
#     'Rn': 10.75, 'Fr': 4.07, 'Ra': 5.28, 'Ac': 5.17, 'Th': 6.31, 
#     'Pa': 5.89, 'U': 6.19, 'Np': 6.27, 'Pu': 6.03, 'Am': 5.97, 
#     'Cm': 5.99, 'Bk': 6.20, 'Cf': 6.28, 'Es': 6.42, 'Fm': 6.50, 
#     'Md': 6.58, 'No': 6.65, 'Lr': 4.90, 'Rf': None, 'Db': None, 
#     'Sg': None, 'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 
#     'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None, 'Mc': None, 
#     'Lv': None, 'Ts': None, 'Og': None
# }

eln = {'D': 13.60, 'OH': 100,
    "H": 2.20, "He": 4.16, "Li": 0.98, "Be": 1.57, "B": 2.04,
    "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98, "Ne": 4.79,
    "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19,
    "S": 2.58, "Cl": 3.16, "Ar": 3.24, "K": 0.82, "Ca": 1.00,
    "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55,
    "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96,
    "Kr": 3.00, "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33,
    "Nb": 1.6, "Mo": 2.16, "Tc": 1.9, "Ru": 2.2, "Rh": 2.28,
    "Pd": 2.20, "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96,
    "Sb": 2.05, "Te": 2.1, "I": 2.66, "Xe": 2.6, "Cs": 0.79,
    "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14,
    "Pm": 1.13, "Sm": 1.17, "Eu": 1.2, "Gd": 1.2, "Tb": 1.2,
    "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Yb": 1.1,
    "Lu": 1.27, "Hf": 1.3, "Ta": 1.5, "W": 2.36, "Re": 1.9,
    "Os": 2.2, "Ir": 2.20, "Pt": 2.28, "Au": 2.54, "Hg": 2.00,
    "Tl": 1.62, "Pb": 2.33, "Bi": 2.02, "Po": 2.0, "At": 2.2,
    "Rn": 2.2, "Fr": 0.7, "Ra": 0.9, "Ac": 1.1, "Th": 1.3,
    "Pa": 1.5, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.3,
    "Cm": 1.3, "Bk": 1.3, "Cf": 1.3, "Es": 1.3, "Fm": 1.3,
    "Md": 1.3, "No": 1.3, "Lr": 1.3, "Rf": 1.3, "Db": 1.3,
    "Sg": 1.3, "Bh": 1.3, "Hs": 1.3, "Mt": 1.3, "Ds": 1.3,
    "Rg": 1.3, "Cn": 1.3, "Nh": 1.3, "Fl": 1.3, "Mc": 1.3,
    "Lv": 1.3, "Ts": 1.3, "Og": 1.3
}

alkali = ['Li', 'Na', 'K', 'Rb']
actinides = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
lanthanides = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
exceptions = actinides + lanthanides

with open(os.path.dirname(__file__)+"/shannon-radii.json") as f:
        out = f.read()
table_ir = json.loads(out)

rn = ["N", "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]