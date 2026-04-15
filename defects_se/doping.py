from pymatgen.analysis.local_env import  CrystalNN
from siman.geo import supercell
import numpy as np 
from .mat_data import  eln, exceptions, table_ir, ionic_radius
from pymatgen.io.vasp import Chgcar
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition
import os
from pymatgen.analysis.defects.generators import ChargeInterstitialGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from sevenn.sevennet_calculator import SevenNetCalculator

from siman.core.structure import Structure
pymatgen2siman = Structure().update_from_pymatgen


class Dopants:
   def __init__(self, alk, st =None, chgcar_file = None, mp_id = None, Int = False):
      """
      NNPUT:
         st - siman structure object
         alk (str) - alkali metal

      OUTPUT:
         class attributes:
            self.alk (str) - alkali metal
            self.dict_oxi (dict) - oxi state by elements
            self.ion_types (dict) - type of ions and their properties  ('coord_number' - coordination, 'non_sym' - non equivalent positions, ...)
      """
      self.alk = alk

      if st:
         self.st = st
         stp = st.convert2pymatgen()

      chgcar = None
      if mp_id:
         try:
            # Access the environment variable named 'API_KEY'
            MP_API_KEY = os.environ['MP_API_KEY']
         except KeyError:
            print("Error: MP_API_KEY environment variable not set.")
         with MPRester(MP_API_KEY) as mpr:
            chgcar = mpr.get_charge_density_from_material_id(mp_id)

               
      if chgcar_file:
         chgcar = Chgcar.from_file(chgcar)

      if chgcar:
         stp =  chgcar.structure
         st = pymatgen2siman(stp)
         self.st = st.copy()
         if Int:
            cig = ChargeInterstitialGenerator(max_insertions = 1)
            ins = cig.generate(chgcar, insert_species=["Pu"])
            defect = next(ins)
            self.sti = self.st.add_atom(defect.site.frac_coords, 'Pu')
         
      at_list =  st.get_unique_type_els()
      oxi_list = st.get_oxi_states('guess')
      
      
      dict_oxi = {}
      for i, el in enumerate(st.get_elements()):
          dict_oxi[el] = oxi_list[i]
      self.dict_oxi = dict_oxi
      
      self.ion_types = {}
      for at in at_list:
          ion_dict = {}
          if dict_oxi[at]>0:
             ion_dict['oxi_state'] = dict_oxi[at]

             sym_pos_ll = st.determine_symmetry_positions(at)
             non_sp = []
             for pos in sym_pos_ll:
                non_sp.append(pos[0])
             ion_dict['non_sym'] =   non_sp

             coord_number = []
             for ind in non_sp:
                neighbors = CrystalNN().get_nn_info(stp, ind)
                coord = len(neighbors)
                coord_number.append(coord)
             ion_dict['coord_number'] = coord_number

             ir_list = []
             for i in range(len(non_sp)):
                 ir = ionic_radius(at, str(int(dict_oxi[at])), str(coord_number[i]))
                 ir_list.append(ir)
             ion_dict['ionic_radius'] = ir_list

             av_dist = []
             for i in range(len(non_sp)):
                av_dist.append(round(st.nn(non_sp[i], coord_number[i])['av'], 3))
             ion_dict['av_dist'] = av_dist

             self.ion_types[at] = ion_dict

            

   def look_match(self):
      """
      INPUT:
        self
   
      UTPUT:
        DR (float) - the smallest relative difference in ionic radii,
        at - atoms for which it is observed
   
      """
      DR = []
      at = []
      for at1 in self.ion_types:
         for at2 in self.ion_types:
            n1 = len(self.ion_types[at1]['non_sym'])
            n2 = len(self.ion_types[at2]['non_sym'])
            for i in range (n1):
               for j in range (n2):
                   if ((at1 != at2) and (self.ion_types[at1]['oxi_state'] != self.ion_types[at2]['oxi_state'])):
                       DR.append(abs(self.ion_types[at1]['ionic_radius'][i] - self.ion_types[at2]['ionic_radius'][j])/ \
                       max(self.ion_types[at1]['ionic_radius'][i], self.ion_types[at2]['ionic_radius'][j]))
                       at.append([at1, at2])
      try:
        ind = DR.index(min(DR))
        self.match_atoms = (min(DR), at[ind])
      except:
         self.match_atoms = None

      return self.match_atoms
       
       

   def search_dop(self, prec_r, prec_eln, d_type, reduce = 1):   
      """
      INPUT:
        prec_r (float) - search dopants with ionic radius in interval  R_i - prec_r * R_i < R_d < R_i + prec_r * R_i
        prec_eln (float) - search dopants with abs(electronegativity diference) < prec_eln
        d_type ( 'vac', 'int') -
             if 'vac' dif = 1 - search dopants with oxi_state(dopant) = oxi_state(ion) + 1
             if 'int' dif = -1 - search dopants with oxi_state(dopant) = oxi_state(ion) - 1
        reduce -
             1 - analyze size of same coordination sphere and serve the best one
             2 - analyze posibility of charge compensation
        OUTPUT:
           list of possible dopants dict
      """
      if d_type == 'vac':
         dif = 1
      elif d_type=='int':
         dif = -1
   
      dop = {}
      for d_el, d_el_dict  in table_ir.items():
         if d_el not in exceptions:
            for d_oxi, d_oxi_dict in d_el_dict.items():
               for d_coord, d_radii in d_oxi_dict.items():
                  d_oxi = float(d_oxi)
                  d_coord = int(d_coord)
                  d_eln = eln[d_el]
                  for ion in self.ion_types:
                        m_eln = eln[ion]
                        n = len(self.ion_types[ion]['non_sym'])
                        for i in range(n):
                             
                           M_IR = self.ion_types[ion]['ionic_radius'][i]
                           r_min = (1-prec_r) * M_IR
                           r_max = (1+prec_r) * M_IR

                           if (
                               (int(self.ion_types[ion]['coord_number'][i]) == d_coord) and 
                               (d_radii>r_min and d_radii<r_max) and
                               (abs(m_eln - d_eln) < prec_eln)  and (d_el not in self.ion_types)  and
                               (self.ion_types[ion]['oxi_state'] + dif == d_oxi)
                              ):

                              key = f'{d_el}_{ion}_{str(i)}'
                              dop_s = {}
                              dop_s['dopant'] = d_el
                              dop_s['matrix_el'] = ion
                              dop_s['matrix_oxi'] = self.ion_types[ion]['oxi_state']
                              dop_s['dop_oxi'] = d_oxi
                              dop_s['matrix_ir'] = M_IR
                              dop_s['dop_ir'] = d_radii
                              dr = (d_radii - M_IR)/M_IR
                              dr = round(dr, 2)
                              dop_s['Delta_R_DM'] = dr
                              dop_s['Delta_ELN_DM'] = round(abs(m_eln - d_eln), 2)
                              dop_s['dop_coord'] = d_coord
                              dop_s['dop_position'] = self.ion_types[ion]['non_sym'][i]
                              dop_s['av_dist'] = self.ion_types[ion]['av_dist'][i]
                              dop[key] = dop_s
                              print(dop)
      dop2 = {}
      if reduce > 0:
         curent_dict = {}
         # orig_key = {}
         for key0, item in dop.items():
            key = (item['matrix_el'], item['dopant'], item['dop_coord'])
            if key not in curent_dict:
               curent_dict[key] = item
               # orig_key[key] = key0
            else:
               if item['Delta_R_DM'] > 0:
                  if item['av_dist'] > curent_dict[key]['av_dist']:
                     curent_dict[key] = item
                     # orig_key[key] = key0
               else:
                  if item['av_dist'] < curent_dict[key]['av_dist']:
                     curent_dict[key] = item
                     # orig_key[key] = key0

         new_dop = {}
         for key, item in curent_dict.items():
             new_dop[f"{item['matrix_el']}_{item['dopant']}_coord{item['dop_coord']}"] = item
         dop = new_dop


      if reduce == 2:
         dop2 = {}
         try:
            rr = self.look_match()
            if rr:
               if rr[0]< 0.05:

                  compare = [rr[1][0], rr[1][1]]
                  if dif == 1:
                     atom = min(compare, key=lambda k: self.ion_types[k]['oxi_state'])
                  if dif == -1:
                     atom = max(compare, key=lambda k: self.ion_types[k]['oxi_state'])

                  for key0, item in dop.items():
                     if item['matrix_el'] != atom:
                        dop2[key0] = item
                  dop = dop2
         except:
            pass
        
      self.dopants = dop
      return(dop)
   

class Dop_Vac:
   def __init__(self, d_obj, key, sc_size):
            

      """
      INPUT:
         d_obj (str) - Dopants object
         key - key of dopant_dict from search_dop
         sc_size - supercell size in A in one direction    
      """
      if max(d_obj.st.rprimd_len()) < sc_size:
         l = [sc_size, sc_size, sc_size]
      else:
         l = []
         for i in range(3):
            if self.st.rprimd_len()[i]<sc_size:
               l.append(sc_size)
            else:
               l.append(d_obj.st.rprimd_len()[i])

      sc = supercell(d_obj.st, l,  test_natom=0)
      self.sc = sc
               

      in_n = self.sc.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
      for j in range(len(in_n)):
         if in_n[j]==d_obj.dopants[key]['dop_position']:
            new_pos = j
            break
      
      dopant = d_obj.dopants[key]['dopant']
      sc_dop = sc.replace_atoms([new_pos], dopant) 
      self.sc_dop = sc_dop

      st_d = []
      r = []
      xred = [] 
      dop_pos =   sc_dop.get_elements_by_el_name(dopant)[0]
      non_sym_alk = sc_dop.determine_symmetry_positions(d_obj.alk)
      nsa = []
      for pos in non_sym_alk:
         nsa.append(pos[0])
      for p in nsa:
         r.append(sc_dop.distance(p, dop_pos))
         xred.append(sc_dop.xred[p])
         cst1 = sc_dop.del_atom(p)
         st_d.append(cst1 )  
      r_st = sorted(zip(r, st_d, xred), key=lambda x: x[0])
      r_sorted, st_sorted, xred_sorted = map(list, zip(*r_st))

      calc = SevenNetCalculator("7net-0") 
      st_short_mine = st_sorted[:3]
      st_dp_mine = [st.convert2pymatgen() for st in st_short_mine]
      st_da_mine = [AseAtomsAdaptor.get_atoms(s) for s in st_dp_mine]
      # sevennet_0_cal = SevenNetCalculator("7net-0")
      E_UP_mine = []
      for s in st_da_mine:
         s.calc = calc
         E_UP_mine.append(s.get_potential_energy())
      ind_min_e = E_UP_mine.index(min(E_UP_mine))
      # st_min_e = st_sorted[ind_min_e]

      self.sc_associate = st_sorted[ind_min_e]
      self.sc_associate_draw = st_sorted[0].add_atom(xred_sorted[ind_min_e], 'Pu')

      self.sc_dissociate = st_sorted[-1]
      self.sc_dissociate_draw = st_sorted[-1].add_atom(xred_sorted[-1], 'Pu')
      self.decompose = Decompose(self.sc_associate, self.sc)



class Dop_Int:
   
   def __init__(self, d_obj, key, sc_size):
            
      """
      INPUT:
         d_obj (str) - Dopants object
         key - key of dopant_dict from search_dop
         sc_size - supercell size in A in one direction    
      """
      if max(d_obj.st.rprimd_len()) < sc_size:
         l = [sc_size, sc_size, sc_size]
      else:
         l = []
         for i in range(3):
            if self.st.rprimd_len()[i]<sc_size:
               l.append(sc_size)
            else:
               l.append(self.st.rprimd_len()[i])

      sc = supercell(d_obj.sti, l, test_natom = 0, skip_numbers = [d_obj.sti.natom - 1])
      self.sc = supercell(d_obj.st, l, test_natom = 0)

      in_n = sc.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
      dop_pos = []
      for j in range(len(in_n)):
         if in_n[j]==d_obj.dopants[key]['dop_position']:
            dop_pos.append(j)

      in_list = sc.get_el('Pu')
      r_min = 100
      r_max = 0
      dpmin = None
      dpmax = None
      for dp in dop_pos:
         r = sc.distance(in_list[0], dp)
         if r>r_max:
            r_max = r
            dpmax = dp
         if r<r_min:
            r_min = r
            dpmin = dp

      dopant = d_obj.dopants[key]['dopant']  
      self.sc_associate_draw = sc.replace_atoms([dpmin], dopant)
      self.sc_dissociate_draw = sc.replace_atoms([dpmax], dopant)
      
      inst_pos = self.sc_associate_draw.get_el('Pu')[0]
      self.sc_associate = self.sc_associate_draw.replace_atoms([inst_pos], d_obj.alk )
      self.sc_dop = self.sc_associate_draw.del_atom(inst_pos)

      inst_pos = self.sc_dissociate_draw.get_el('Pu')[0]
      self.sc_dissociate = self.sc_dissociate_draw.replace_atoms([inst_pos], d_obj.alk )

      self.decompose = Decompose(self.sc_associate, self.sc)      



class Decompose:

    def __init__(self, st, sc):

        MP_API_KEY = os.environ["MP_API_KEY"]

        self.init_st = st

        # химическая система
        ch = "-".join(st.get_unique_type_els())

        with MPRester(MP_API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(
                ch,
                compatible_only=False,
                additional_criteria={
                    "is_stable": True,
                    "thermo_types": ["GGA_GGA+U"],
                },
            )

        # phase diagram
        phd = PhaseDiagram(entries)

        comp_def = Composition(st.get_formula())
        dec_def = phd.get_decomposition(comp_def)

        # список элементов
        elements = sorted(set(st.get_unique_type_els()))

        # функция атомного вектора
        def atom_vector(struct):
            comp = Composition(struct.get_formula())
            return np.array([comp[el] for el in elements])

        # вектор defect
        b = atom_vector(st)

        # список фаз
        phases = []
        phase_vecs = []

        for entry in dec_def:

            st_d = pymatgen2siman(entry.structure)

            # sc будет отдельной переменной
            if st_d.get_reduced_formula() == sc.get_reduced_formula():
                continue

            phases.append(entry)
            phase_vecs.append(atom_vector(st_d))

        # матрица системы
        A = np.column_stack([atom_vector(sc)] + phase_vecs)

        # решение линейной системы
        x, *_ = np.linalg.lstsq(A, b, rcond=None)

        alpha = x[0]
        beta = x[1:]

        # результаты
        self.dec_ids = []
        self.dec_st = []
        self.dec_coef = []

        reaction = st.get_formula() + " -> "

        # host
        reaction += f"{round(alpha,3)} {sc.get_formula()} + "

        self.dec_ids.append("sc")
        self.dec_st.append(sc.copy())
        self.dec_coef.append(alpha)

        # остальные фазы
        for entry, c in zip(phases, beta):

            if abs(c) < 1e-6:
                continue

            st_d = pymatgen2siman(entry.structure)

            self.dec_ids.append(entry.data["material_id"])
            self.dec_st.append(st_d)
            self.dec_coef.append(c)

            reaction += f"{round(c,3)} {entry.formula} + "

        self.dec_reaction = reaction[:-3]



def Ch(st, E0, Eb, matrix, T = 1000, ALK = 'Li', g = 6):
    """
    INPUT:
       st - siman structure
       E0 - decomposition energy(eV)
       Eb - binding energy(eV)
       matrix - matrix element
    OUTPUT:
       conc (float) - concentration of dopants (N_dopants/N_matrix)
    """
    k = 8.617e-05
    NA = len(st.get_elements_by_el_name(ALK))
    M = len(st.get_elements_by_el_name(matrix))
    ml = M/NA

    C = g * np.exp(-E0/(k*T)) + np.exp(-(E0+Eb)/(2*k*T))/np.sqrt(ml)
    return C

def Cl(st, ch, Eb, T=300, g = 6):
    """
    INPUT:
       ch - hight temperature dopant concentration (c_vac = N_matrix/N_alk * c_dop)
       Eb - binding energy(eV)
    OUTPUT:
       conc (float) - free RT vacancy concentration 
    """
    NA = len(st.get_elements_by_el_name(ALK))
    M = len(st.get_elements_by_el_name(matrix))
    ml = M/NA 
    c = ml * ch
    k = 8.617e-05
    al = np.exp(-Eb/(k*T))/g
    C = (np.sqrt(al**2 +4*al*ch) - al) / (2 * (1 - al))
    return C

