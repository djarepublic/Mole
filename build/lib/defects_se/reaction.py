import numpy as np

API_KEY = 'AaqKFXh9Qzoo1FqCIbELtEVHa72BW3EG'
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import (CompoundPhaseDiagram, PDPlotter, PhaseDiagram, PDEntry)
from pymatgen.core import Composition
from siman.core.structure import Structure()
pymatgen2siman = Structure().update_from_pymatgen

class Decompose():
    def __init__(self, st):
      """
      INPUT:
        structure

      OUTPUT:
        list of decomposition mp_id
      """
      self.init_st = st
      ch = st.get_unique_type_els()[0]
      for el in st.get_unique_type_els()[1:]:
          ch = ch + '-' + el    

      with MPRester(API_KEY) as mpr:
           entries = mpr.get_entries_in_chemsys(ch,
                                                  compatible_only=False,
                                                  additional_criteria={'is_stable': True,
                                                                      'thermo_types': ['GGA_GGA+U']})

      phd = PhaseDiagram(entries)
      cst = Composition(st.get_formula())
      dec = phd.get_decomposition(cst)
      ids = []
      dec_st = []
      coef = []
      prod = cst.formula + ' -> '
      for key, val  in dec.items():
          st_d = pymatgen2siman(key.structure)
          st_d.name = key.data['material_id'] + '_' + key.formula
          c = val*st.natom/st_d.natom
          prod = prod +  rond(c, 3) + ' ' + key.formula + ' + '
          ids.append(key.data['material_id'])
          dec_st.append(st_d)
          coef.append(c)
      self.dec_ids = ids
      self.dec_st = dec_st
      self.dec_cof = coef
      self.dec_reaction = prod[:-3]


def find_decomp_id(st):
      """
      INPUT:
        structure

      OUTPUT:
        list of decomposition mp_id
      """
      ch = st.get_unique_type_els()[0]
      for el in st.get_unique_type_els()[1:]:
          ch = ch + '-' + el    

      with MPRester(API_KEY) as mpr:
           entries = mpr.get_entries_in_chemsys(ch,
                                                  compatible_only=False,
                                                  additional_criteria={'is_stable': True,
                                                                      'thermo_types': ['GGA_GGA+U']})

      phd = PhaseDiagram(entries)
      cst = Composition(st.get_formula())
      dec = phd.get_decomposition(cst)
      id_val = []
      for key, val  in dec.items():
          id_val.append((key.data['material_id'], val))

      return id_val


def calc_react(initial_phase, prod):
    """
    INPUT:
        initial_phase - tuple - (initial_structure, energy)
        prod - list of product tuple [(product1, energy1), (product2, energy2), ...]
    OUTPUT:
        str of reaction, energy (per cell)
    """
    init_st = initial_phase[0]
    atoms = init_st.get_unique_type_els()
    matrix = np.zeros([len(prod), len(prod)])
    nubers_inr = []
    for i in range(len(prod)):
        nubers_inr.append(len(init_st.get_elements_by_el_name(atoms[i])))
        numbers_inp = []
        for pr in prod:
            n = len(pr[0].get_elements_by_el_name(atoms[i]))
            numbers_inp.append(n)
        numbers_inp = np.array(numbers_inp)
        matrix[i] = numbers_inp

    nubers_inr = np.array(nubers_inr)
    coef =  np.linalg.solve(matrix, nubers_inr)

    E = initial_phase[1]
    for i in range(len(prod)):
        E = E - coef[i]*prod[i][1]

    rection = initial_phase[0].get_formula() + ' -> '
    for i in range(len(prod)):
        rection = rection + str(coef[i]) + ' ' + prod[i][0].get_formula() + ' + '
    rection = rection[:-3]

    return(rection, E)

        
        
