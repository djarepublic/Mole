
import siman 
from siman.calc_manage import smart_structure_read, get_structure_from_matproj,  res_loop
import json    
from siman.mat_prop.mat_prop import ionic_radius
from pymatgen.analysis.local_env import VoronoiNN, CrystalNN
from siman.geo import supercell
from pymatgen.analysis.energy_models import EwaldElectrostaticModel

from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from mp_api.client import MPRester
from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import (CompoundPhaseDiagram, PDPlotter, PhaseDiagram, PDEntry)
from siman.core.structure import Structure
import numpy as np 

rn = ["N", "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]

with open("/home/dany/research/skoltech/doping/mat_prop/shannon-radii.json") as f:
        out = f.read()
table_ir = json.loads(out)

ionization_energies = { 'D': 13.60, 'OH': 100,
    'H': 13.60, 'He': 24.59, 'Li': 5.39, 'Be': 9.32, 'B': 8.30, 
    'C': 11.26, 'N': 14.53, 'O': 13.62, 'F': 17.42, 'Ne': 21.56, 
    'Na': 5.14, 'Mg': 7.65, 'Al': 5.99, 'Si': 8.15, 'P': 10.49, 
    'S': 10.36, 'Cl': 12.97, 'Ar': 15.76, 'K': 4.34, 'Ca': 6.11, 
    'Sc': 6.56, 'Ti': 6.83, 'V': 6.75, 'Cr': 6.77, 'Mn': 7.43, 
    'Fe': 7.90, 'Co': 7.88, 'Ni': 7.64, 'Cu': 7.73, 'Zn': 9.39, 
    'Ga': 6.00, 'Ge': 7.90, 'As': 9.79, 'Se': 9.75, 'Br': 11.81, 
    'Kr': 14.00, 'Rb': 4.18, 'Sr': 5.69, 'Y': 6.22, 'Zr': 6.63, 
    'Nb': 6.76, 'Mo': 7.09, 'Tc': 7.28, 'Ru': 7.36, 'Rh': 7.46, 
    'Pd': 8.34, 'Ag': 7.58, 'Cd': 8.99, 'In': 5.79, 'Sn': 7.34, 
    'Sb': 8.61, 'Te': 9.01, 'I': 10.45, 'Xe': 12.13, 'Cs': 3.89, 
    'Ba': 5.21, 'La': 5.58, 'Ce': 5.54, 'Pr': 5.47, 'Nd': 5.53, 
    'Pm': 5.58, 'Sm': 5.64, 'Eu': 5.67, 'Gd': 6.15, 'Tb': 5.86, 
    'Dy': 5.94, 'Ho': 6.02, 'Er': 6.11, 'Tm': 6.18, 'Yb': 6.25, 
    'Lu': 5.43, 'Hf': 6.83, 'Ta': 7.55, 'W': 7.86, 'Re': 7.83, 
    'Os': 8.44, 'Ir': 8.97, 'Pt': 8.96, 'Au': 9.23, 'Hg': 10.44, 
    'Tl': 6.11, 'Pb': 7.42, 'Bi': 7.29, 'Po': 8.42, 'At': 9.30, 
    'Rn': 10.75, 'Fr': 4.07, 'Ra': 5.28, 'Ac': 5.17, 'Th': 6.31, 
    'Pa': 5.89, 'U': 6.19, 'Np': 6.27, 'Pu': 6.03, 'Am': 5.97, 
    'Cm': 5.99, 'Bk': 6.20, 'Cf': 6.28, 'Es': 6.42, 'Fm': 6.50, 
    'Md': 6.58, 'No': 6.65, 'Lr': 4.90, 'Rf': None, 'Db': None, 
    'Sg': None, 'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 
    'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None, 'Mc': None, 
    'Lv': None, 'Ts': None, 'Og': None
}

alkali = ['Li', 'Na', 'K', 'Rb']
actinides = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
lanthanides = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
ex = actinides + lanthanides

class defects:
    def __init__(self, st, alk):
        """
      NNPUT:
         st - siman structure object
         alk (str) - alkali metal
      
      RETURN:
         class attributes:
            self.alk (str) - alkali metal
            self.dict_oxi (dict) - oxi state by elements
            self.ion_types (dict) - type of ions and their properties  ('ip' - ionizational potential, 'non_sym' - non equivalent positions, ...)
        """
        self.alk = alk
        at_list = st.get_unique_type_els()
        oxi_list = st.get_oxi_states('guess')
        stp = st.convert2pymatgen()
        stp.add_oxidation_state_by_guess()
        
        dict_oxi = {}
        for i, el in enumerate(st.get_elements()):
            dict_oxi[el] = oxi_list[i]
        self.dict_oxi = dict_oxi

        self.st = st
        self.ion_types = {}

        for at in at_list:
            #self.ion_types[at] = {}
            ion_dict = {}
            ion_dict['ip'] =  ionization_energies[at] #in eV
    
            sym_pos_ll = st.determine_symmetry_positions(at)
            l_sp = []
            for pos in sym_pos_ll:
               l_sp.append(pos[0])
            ion_dict['non_sym'] =   l_sp
            #ion_dict['number_pos'] = len(l_sp)
            oxi_list_at = []
            for ind in l_sp:
               oxi_list_at.append(oxi_list[ind])
            ion_dict['oxi_state'] = oxi_list_at
            
            coord_number = []
            for ind in l_sp:
               neighbors = CrystalNN().get_nn_info(stp, ind)
               coord = len(neighbors)
               coord_number.append(coord)
            ion_dict['coord_number'] = coord_number
    
            ir_list = []
            for i in range(len(l_sp)):
               
               try:
                ir = ionic_radius(at, int(oxi_list_at[i]), rn[coord_number[i]])
               except:
                for j in range(4):
                  try:
                     ir = ionic_radius(at, int(oxi_list_at[i]), rn[coord_number[i]+j])
                  except:
                     continue
                  if ir:
                     break
               ir_list.append(ir)
            ion_dict['ionic_radius'] = ir_list

            if ion_dict['oxi_state'][0] > 0:
               self.ion_types[at] = ion_dict

            

    def look_match(self, prec_rr=0.03):
       """
       INPUT:
         self
      
      OUTPUT:
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
                    if ((at1 != at2) and (self.ion_types[at1]['oxi_state'][i] != self.ion_types[at2]['oxi_state'][j])):
                        DR.append(abs(self.ion_types[at1]['ionic_radius'][i] - self.ion_types[at2]['ionic_radius'][j])/ \
                        max(self.ion_types[at1]['ionic_radius'][i], self.ion_types[at2]['ionic_radius'][j]))
                        at.append([at1, at2])
       ind = DR.index(min(DR))
         #                a = 1
         #                break
         #        if a==1:
         #           break
         #     if a==1:
         #      break
         #  if a==1:
         #     break
       return (min(DR), at[ind])
       
       

    def search_dop(self, prec_r, prec_ip, dif=None):
       """
       INPUT:
         prec_r (float) - search dopants with ionic radius in interval  R_i - prec_r * R_i < R_d < R_i + prec_r * R_i
         prec_ip (float) - search dopants with ionization potential in interval  IP_i - prec_ip * IP_i < IP_d < IP_i + prec_ip * IP_i
         dif (None, 1, 2) -
              if 1 - search dopants with oxi_state(dopant) = oxi_state(ion) + 1
              if 2 - search dopants with oxi_state(dopant) = oxi_state(ion) + 2
              if None - both cases
         OUTPUT:
            dopants dict
       """
       dop = {}
       for d_el, d_el_dict  in table_ir.items():
         dop_s = {}
         d_oxi_list = []
         d_coord_list = []
         d_ir_list = []
         d_rep_elements = []
         d_position = []
         el_oxi_list = []
         for d_oxi, d_oxi_dict in d_el_dict.items():
             for d_coord, prop in d_oxi_dict.items():
                 if (d_coord == "IVSQ" or d_coord == "IVPY"):
                     d_coord = "IV"
                 if (d_coord == "VILS" or d_coord == "VIHS"):
                     d_coord = "VI"
                 if (d_coord == 'IIIPY'):
                     d_coord = 'III'
                 d_oxi = float(d_oxi)
                 d_coord = int(rn.index(d_coord))
                 d_ip = ionization_energies[d_el]
                 for ion in self.ion_types:
                        n = len(self.ion_types[ion]['non_sym'])
                        for i in range(n):
                           try:
                            r_min = (1-prec_r)*self.ion_types[ion]['ionic_radius'][i]
                            r_max = (1+prec_r)*self.ion_types[ion]['ionic_radius'][i]
                            ip_min = (1-prec_ip)*self.ion_types[ion]['ip']
                            ip_max = (1+prec_ip)*self.ion_types[ion]['ip']
                      
                            if (
                                (int(self.ion_types[ion]['coord_number'][i]) == d_coord) and 
                                (prop['r_ionic']>r_min and prop['r_ionic']<r_max) and
                                (d_el not in dop) and (d_ip > ip_min and d_ip < ip_max) and (d_el not in ex)
                               ):
                              
                              if( 
                                 ((dif==1) and (self.ion_types[ion]['oxi_state'][i] +1 == d_oxi )) or 
                                  ((dif==2) and (self.ion_types[ion]['oxi_state'][i] +2 == d_oxi )) or
                                  ((self.ion_types[ion]['oxi_state'][i] +1 == d_oxi or self.ion_types[ion]['oxi_state'][i] +2 == d_oxi) and (dif==None))
                                 ):
                               
                               el_oxi_list.append(self.ion_types[ion]['oxi_state'][i])
                               d_oxi_list.append(d_oxi)
                               d_coord_list.append(d_coord)
                               d_rep_elements.append(ion)
                               d_ir_list.append(prop['r_ionic'])
                               d_position.append(self.ion_types[ion]['non_sym'][i])
                           
                           except TypeError as e:
                              continue
                           

         if d_coord_list:
            dop_s['elements'] = d_rep_elements
            dop_s['position'] = d_position
            dop_s['oxi_state'] = d_oxi_list
            dop_s['el_oxi_states'] = el_oxi_list
            dop_s['coord_number'] = d_coord_list
            dop_s['ir'] = d_ir_list
         if dop_s:
            dop[d_el] = dop_s
       return(dop)
            


    def make_st(self, dop_el, dop_dict, C = 8, st_numb = 10):
        """
         INPUT:
            dop_el (str) - dopant element
            dop_dict (dict) - dopant_dict from search_dop
            C (int) - reverse concentration
            st_numb (int) - number of output structures
         OUTPUT:
            st (list) - list of structure sorted by increasing Ewald energy
            ew (list) - list of Ewald energies sorted in ascending order
            r (list) - list if radii btw dopant and vacancy sorted by increasing Ewald energy (only if one vacancy)

        """
        
        st_d = []
        
        for i in range(len(dop_dict['elements'])):
            NN = C
            if dop_dict['elements'][i] == self.alk:
               if dop_dict['oxi_state'][i] - dop_dict['el_oxi_states'][i] == 2:
                  NN = NN + 3
               else:
                  NN = NN + 2
            # else:
            #    if dop_dict['oxi_state'][i] - dop_dict['el_oxi_states'][i] == 2:
            #       NN = NN + 2
            #    else:
            #       NN = NN + 1 
            N = len(self.st.get_numbers(dop_dict['elements'][i]))
            l = 10
            c_st = supercell(self.st, [l, l, l])
            
            while N < NN:
               l = l + 1
               c_st = supercell(self.st, [l, l, l])
               N = len(c_st.get_numbers(dop_dict['elements'][i]))
            # print(l)
            # k = 1.380649 * 10**(-23)
            # conf_entropy = k*np.log(1/N)
            # T = 1000
            in_n = c_st.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
            for j in range(len(in_n)):
               if in_n[j]==dop_dict['position'][i]:
                  new_pos = j
                  break

            
            

            r = []
            if dop_dict['oxi_state'][i] - dop_dict['el_oxi_states'][i] == 1:
               cst = c_st.replace_atoms([new_pos], dop_el)    
               dop_pos =   cst.get_elements_by_el_name(dop_el)[0]
               non_sym_alk = cst.determine_symmetry_positions(self.alk)
               nsa = []
               for pos in non_sym_alk:
                  nsa.append(pos[0])
               for p in nsa:
                  r.append(cst.distance(p, dop_pos))
                  cst1 = cst.del_atom(p)
                  st_d.append(cst1)

            
            if dop_dict['oxi_state'][i] - dop_dict['el_oxi_states'][i] == 2:
               cst = c_st.replace_atoms([new_pos], dop_el)
               non_sym_alk = cst.determine_symmetry_positions(self.alk)
               nsa = []
               for pos in non_sym_alk:
                  nsa.append(pos[0])
               for p in nsa:
                  cst1 = cst.del_atom(p)
                  non_sym_alk2 = cst1.determine_symmetry_positions(self.alk)
                  nsa2 = []
                  for pos in non_sym_alk2:
                     nsa2.append(pos[0])
                  for p2 in nsa2:
                     cst2 = cst1.del_atom(p2)
                     if cst2 not in st_d:
                        st_d.append(cst2)

        st_dp_unic = [st.convert2pymatgen() for st in st_d]
   #      st_dp_unic = []
   #      for structure in st_dp:
   #  # Проверяем, есть ли уже такая структура в unique_structures
   #          if not any(structure == existing for existing in st_dp_unic):
   #             st_dp_unic.append(structure)

      #   if ew==False:
      #    return st_d
      #        else:
        st_e = []            
        for stp in st_dp_unic:
           dict_oxi = self.dict_oxi
           dict_oxi[dop_el] = dop_dict['oxi_state'][0]
           stp.add_oxidation_state_by_element(dict_oxi)
           E_ev = EwaldElectrostaticModel().get_energy(stp)
           st = Structure().update_from_pymatgen(stp)
           st_e.append((st, E_ev))      
        st_e = sorted(st_e, key=lambda list: list[1])
        st = []
        ew = []
        st_numb = min(st_numb, len(st_e))
        for i in range(st_numb):
         st.append(st_e[i][0])
         ew.append(st_e[i][1])
        ew = [e - min(ew) for e in ew]
        
        return st, ew, r


    def decomp(self, calc, API_KEY = 'AaqKFXh9Qzoo1FqCIbELtEVHa72BW3EG'):
      """
      INPUT:
         calc - siman calc
      OUTPUT:
         react_E (float) - energy of reaction above hull
         reaction_C (tuple) - reaction as (str) and coefficints
         conc (float) - concentration 

      """

      

      ch = calc.end.get_unique_type_els()[0]
      for el in calc.end.get_unique_type_els()[1:]:
        ch = ch + '-' + el    

      comp = Composition(calc.end.get_formula())
      my_entry = PDEntry(comp, calc.e0)

      with MPRester(API_KEY) as mpr:
         entries = mpr.get_entries_in_chemsys(ch,
                                                    compatible_only=False,
                                                    additional_criteria={'is_stable': True,
                                                                        'thermo_types': ['GGA_GGA+U']})
      reaction = str(calc.end.get_formula()) + ' → '
      pd = PhaseDiagram(entries)
      entries1 = pd.get_decomp_and_phase_separation_energy(my_entry, stable_only = True)

      i = 0 
      C = []
      for key in entries1[0].keys():
         i = i + 1
         reaction = reaction + 'С_' + str(i) +' '+ str(key.composition) + ' + '
         C.append(entries1[0][key])
      reaction = reaction[:-3]

      reaction_C = (reaction, C)
      react_E = entries1[1] * calc.end.natom

      k = 8.6*10**(-5)
      T = 1000 
      conc = np.exp(-react_E/(2*k*T))
      
      return(react_E, reaction_C, conc)


