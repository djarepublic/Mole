import numpy as np
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import (CompoundPhaseDiagram, PDPlotter, PhaseDiagram, PDEntry)
from pymatgen.core import Composition
from siman.core.structure import Structure
pymatgen2siman = Structure().update_from_pymatgen
import os


class Decompose():
    def __init__(self, st):
      """
      INPUT:
        structure

      OUTPUT:
        list of decomposition mp_id
      """
      try:
          # Access the environment variable named 'API_KEY'
          MP_API_KEY = os.environ['MP_API_KEY']
      except KeyError:
          print("Error: MP_API_KEY environment variable not set.")

      self.init_st = st
      ch = st.get_unique_type_els()[0]
      for el in st.get_unique_type_els()[1:]:
          ch = ch + '-' + el    

      with MPRester(MP_API_KEY) as mpr:
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
          prod = prod +  str(round(c, 3)) + ' ' + key.formula + ' + '
          ids.append(key.data['material_id'])
          dec_st.append(st_d)
          coef.append(c)
      self.dec_ids = ids
      self.dec_st = dec_st
      self.dec_cof = coef
      self.dec_reaction = prod[:-3]




        
        
