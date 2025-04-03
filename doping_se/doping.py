
import siman 
from siman.calc_manage import smart_structure_read, get_structure_from_matproj,  res_loop
import json    
from pymatgen.analysis.local_env import VoronoiNN, CrystalNN
from siman.geo import supercell
from pymatgen.analysis.energy_models import EwaldElectrostaticModel

from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from mp_api.client import MPRester
from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import (CompoundPhaseDiagram, PDPlotter, PhaseDiagram, PDEntry)
from siman.core.structure import Structure
import numpy as np 
#from sevenn.sevennet_calculator import SevenNetCalculator
from sevenn.calculator import SevenNetCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from collections import defaultdict

from mat_data import  eln, exceptions, table_ir, rn, ionic_radius
from reaction import find_decomp_id, calc_react


class defects:
    def __init__(self, st, alk):
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
        st = st.get_primitive_cell()
        at_list = st.get_unique_type_els()
        oxi_list = st.get_oxi_states('guess')
        stp = st.convert2pymatgen()
      #   stp.add_oxidation_state_by_guess()
        
        dict_oxi = {}
        for i, el in enumerate(st.get_elements()):
            dict_oxi[el] = oxi_list[i]
        self.dict_oxi = dict_oxi

        self.st = st
        self.ion_types = {}

        for at in at_list:
            #self.ion_types[at] = {}
            ion_dict = {}
            ion_dict['oxi_state'] = dict_oxi[at]
    
            sym_pos_ll = st.determine_symmetry_positions(at)
            l_sp = []
            for pos in sym_pos_ll:
               l_sp.append(pos[0])
            ion_dict['non_sym'] =   l_sp
            
            coord_number = []
            for ind in l_sp:
               neighbors = CrystalNN().get_nn_info(stp, ind)
               coord = len(neighbors)
               coord_number.append(coord)
            ion_dict['coord_number'] = coord_number
    
            ir_list = []
            for i in range(len(l_sp)):
               
               try:
                ir = ionic_radius(at, int(dict_oxi[at]), rn[coord_number[i]])
               except:
                for j in range(4):
                  try:
                     ir = ionic_radius(at, int(dict_oxi[at]), rn[coord_number[i]+j])
                  except:
                     continue
                  if ir:
                     break
               ir_list.append(ir)
            ion_dict['ionic_radius'] = ir_list

            av_dist = []
            for i in range(len(l_sp)):
               av_dist.append(round(st.nn(l_sp[i], coord_number[i])['av'], 3))
            ion_dict['av_dist'] = av_dist

            if ion_dict['oxi_state'] > 0:
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
                    if ((at1 != at2) and (self.ion_types[at1]['oxi_state'] != self.ion_types[at2]['oxi_state'])):
                        DR.append(abs(self.ion_types[at1]['ionic_radius'][i] - self.ion_types[at2]['ionic_radius'][j])/ \
                        max(self.ion_types[at1]['ionic_radius'][i], self.ion_types[at2]['ionic_radius'][j]))
                        at.append([at1, at2])
       ind = DR.index(min(DR))

       return (min(DR), at[ind])
       
       

    def search_dop(self, prec_r, prec_eln, dif=None, reduce = 1):   #prec_ip
       """
       INPUT:
         prec_r (float) - search dopants with ionic radius in interval  R_i - prec_r * R_i < R_d < R_i + prec_r * R_i
         prec_eln (float) - search dopants with abs(electronegativity diference) < prec_eln
         dif (None, 1, 2) -
              if 1 - search dopants with oxi_state(dopant) = oxi_state(ion) + 1
              if 2 - search dopants with oxi_state(dopant) = oxi_state(ion) + 2
              if None - both cases
         reduce -
              1 - analyze size of coordination sphere and serve the best one
              0 - all cases
         OUTPUT:
            list of possible dopants dict
       """
       dop = []
       for d_el, d_el_dict  in table_ir.items():
         # dop_s = {}
         # d_oxi_list = []
         # d_coord_list = []
         # d_ir_list = []
         # d_rep_elements = []
         # d_position = []
         # el_oxi_list = []
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
                 #d_ip = ionization_energies[d_el]
                 d_eln = eln[d_el]
                 for ion in self.ion_types:

                        m_eln = eln[ion]
                        n = len(self.ion_types[ion]['non_sym'])
                        for i in range(n):
                            
                           # try:
                            M_IR = self.ion_types[ion]['ionic_radius'][i]
                            r_min = (1-prec_r) * M_IR
                            r_max = (1+prec_r) * M_IR
                            #ip_min = (1-prec_ip)*self.ion_types[ion]['ip']
                            #ip_max = (1+prec_ip)*self.ion_types[ion]['ip']
                      
                            if (
                                (int(self.ion_types[ion]['coord_number'][i]) == d_coord) and 
                                (prop['r_ionic']>r_min and prop['r_ionic']<r_max) and
                                (abs(m_eln - d_eln) < prec_eln) and (d_el not in exceptions) and
                                (d_el not in self.ion_types) 
                                #(d_ip > ip_min and d_ip < ip_max) and (d_el not in exceptions)
                               ):
                              
                              if( 
                                 ((dif==1) and (self.ion_types[ion]['oxi_state'] +1 == d_oxi )) or 
                                  ((dif==2) and (self.ion_types[ion]['oxi_state'] +2 == d_oxi )) or
                                  ((self.ion_types[ion]['oxi_state'] +1 == d_oxi or self.ion_types[ion]['oxi_state'] +2 == d_oxi) and (dif==None))
                                 ):
                               

                              #  el_oxi_list.append(self.ion_types[ion]['oxi_state'][i])
                              #  d_oxi_list.append(d_oxi)
                               dop_s = {}
                               dop_s['dopant'] = d_el
                               dop_s['matrix_el'] = ion
                               dop_s['matrix_oxi'] = self.ion_types[ion]['oxi_state']
                               dop_s['dop_oxi'] = d_oxi
                               dop_s['matrix_ir'] = M_IR
                               dop_s['dop_ir'] = prop['r_ionic']
                               dr = (prop['r_ionic'] - M_IR)/M_IR
                               dr = round(dr, 2)
                               dop_s['Delta_R_DM'] = dr
                               dop_s['Delta_ELN_DM'] = round(abs(m_eln - d_eln), 2)
                              #  DR.append(dr) 
                               dop_s['dop_coord'] = d_coord
                              #  d_rep_elements.append(ion)
                               
                               dop_s['dop_position'] = self.ion_types[ion]['non_sym'][i]
                               dop_s['av_dist'] = self.ion_types[ion]['av_dist'][i]
                               dop.append(dop_s)
                           
                           # except TypeError as e:
                           #    continue
                           
  
                              #  if d_coord_list:
                              #     dop_s['element'] = rep_el
                              #     dop_s['oxi_state'] = d_oxi
                              #     dop_s['position'] = d_position
                              #     dop_s['DR'] = DR
                              #     # dop_s['el_oxi_states'] = el_oxi_list
                              #     dop_s['coord_number'] = d_coord_list
                              #     dop_s['ir'] = d_ir_list
                              #  if dop_s:
                              #     dop[d_el] = dop_s
       if reduce == 1:
         result_dict = {}

         for item in dop:
            key = (item['matrix_el'], item['dopant'], item['dop_coord'])

            if key not in result_dict:
               result_dict[key] = item
            else:
               if item['Delta_R_DM'] > 0:
                  if item['av_dist'] > result_dict[key]['av_dist']:
                     result_dict[key] = item
               else:
                  if item['av_dist'] < result_dict[key]['av_dist']:
                     result_dict[key] = item
         
         dop = list(result_dict.values())

         dop2 = []
         try:
            rr = self.look_match()
            if rr:
               if rr[0]< 0.1:
                  target_key = 'oxi_state'
                  compare = [rr[1][0], rr[1][1]]
                  #print(compare)
                  atom = min(compare, key=lambda k: self.ion_types[k][target_key])
                  print(atom)
                  for item in dop:
                     if item['matrix_el'] != atom:
                        dop2.append(item)
                  dop = dop2
         except:
            pass
         
       
       return(dop)
            


    def make_st(self, dop_el, dop_dict):
        """
         INPUT:
            dop_el (str) - dopant element
            dop_dict (dict) - dopant_dict from search_dop

         OUTPUT:
            r_sorted - distance between vacancy and dopant
            st_sorted - structures
            xred_sorted - xred of vacancies

        """
        
        st_d = []
        

      #   NN = C
      #   if dop_dict['matrix_el'] == self.alk:
      #      if dop_dict['dop_oxi'] -  dop_dict['matrix_oxi'] == 2:
      #         NN = NN + 2
      #      else:
      #         NN = NN + 1
        # else:
        #    if dop_dict['oxi_state'][i] - dop_dict['el_oxi_states'][i] == 2:
        #       NN = NN + 2
        #    else:
        #       NN = NN + 1 
      #   N = len(self.st.get_numbers(dop_dict['matrix_el']))
        if max(self.st.rprimd_len()) < 10:
            l = [10, 10, 10]
            c_st = supercell(self.st, l)
        else:
            l = []
            for i in range(3):
               if self.st.rprimd_len()[i]<10:
                  l.append(10)
               else:
                  l.append(self.st.rprimd_len()[i])
            c_st = supercell(self.st, l,  test_natom=0)
               



        
      #   while N < NN:
      #      l = [l1 + 1 for l1 in l]
      #      c_st = supercell(self.st, l)
      #      N = len(c_st.get_numbers(dop_dict['matrix_el']))

        in_n = c_st.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
        for j in range(len(in_n)):
           if in_n[j]==dop_dict['dop_position']:
              new_pos = j
              break

            
            

        r = []
        xred = []
        if dop_dict['dop_oxi'] -  dop_dict['matrix_oxi'] == 1:
           cst = c_st.replace_atoms([new_pos], dop_el)    
           dop_pos =   cst.get_elements_by_el_name(dop_el)[0]
           non_sym_alk = cst.determine_symmetry_positions(self.alk)
           nsa = []
           for pos in non_sym_alk:
              nsa.append(pos[0])
           for p in nsa:
              r.append(cst.distance(p, dop_pos))
              xred.append(cst.xred[p])
              cst1 = cst.del_atom(p)
              st_d.append(cst1 ) 

           r_st = sorted(zip(r, st_d, xred), key=lambda x: x[0])
           r_sorted, st_sorted, xred_sorted = map(list, zip(*r_st))
           return r_sorted, st_sorted, xred_sorted


    def make_st2(self, dop_el, dop_dict, draw = False):
            """
             INPUT:
                dop_el (str) - dopant element
                dop_dict (dict) - dopant_dict from search_dop
                draw - add Pu atom in vacancy
    
             OUTPUT:
                st_min_e, st_max_r

                if draw == True:
                  st_min_e, st_max_r, st_min_e + Pu, st_max_r + Pu
    
            """
            
           
            r_sorted, st_sorted, xred_sorted = self.make_st(dop_el, dop_dict)
            st_short_mine = st_sorted[:4]
            ind_max = len(st_sorted) - 1
         #   st_short_maxr = st_sorted[ind_max - 2:]
        
      #   if dop_dict['dop_oxi'] -  dop_dict['matrix_oxi'] == 2:
      #      cst = c_st.replace_atoms([new_pos], dop_el)
      #      non_sym_alk = cst.determine_symmetry_positions(self.alk)
      #      nsa = []
      #      for pos in non_sym_alk:
      #         nsa.append(pos[0])
      #      for p in nsa:
      #         cst1 = cst.del_atom(p)
      #         non_sym_alk2 = cst1.determine_symmetry_positions(self.alk)
      #         nsa2 = []
      #         for pos in non_sym_alk2:
      #            nsa2.append(pos[0])
      #         for p2 in nsa2:
      #            cst2 = cst1.del_atom(p2)
      #            if cst2 not in st_d:
      #               st_d.append(cst2)

            st_dp_mine = [st.convert2pymatgen() for st in st_short_mine]
            st_da_mine = [AseAtomsAdaptor.get_atoms(s) for s in st_dp_mine]

      #       st_dp_maxr = [st.convert2pymatgen() for st in st_short_maxr]
      #       st_da_maxr = [AseAtomsAdaptor.get_atoms(s) for s in st_dp_maxr]


            sevennet_0_cal = SevenNetCalculator("7net-0") 

            E_UP_mine = []
            for s in st_da_mine:
               s.calc = sevennet_0_cal
               E_UP_mine.append(s.get_potential_energy())

      #   E_UP_maxr = []
      #   for s in st_da_maxr:
      #      s.calc = sevennet_0_cal
      #      E_UP_maxr.append(s.get_potential_energy())
        
        
            ind_min_e = E_UP_mine.index(min(E_UP_mine))
      #       ind_max_r = E_UP_maxr.index(min(E_UP_maxr))


            st_min_e = st_sorted[ind_min_e]
            st_max_r = st_sorted[ind_max ]




            if draw:
               st_min_ed = st_min_e.add_atom(xred_sorted[ind_min_e], 'Pu')
               st_max_rd = st_max_r.add_atom(xred_sorted[ind_max ], 'Pu')

               return st_min_e, st_max_r, st_min_ed, st_max_rd

            else:
            
             return st_min_e, st_max_r
           
         
   #      st_dp_unic = []
   #      for structure in st_dp:
   #  # Проверяем, есть ли уже такая структура в unique_structures
   #          if not any(structure == existing for existing in st_dp_unic):
   #             st_dp_unic.append(structure)

      #   if ew==False:
      #    return st_d
      #        else:


      #   ew = []            
      #   for stp in st_dp_unic:
      #      dict_oxi = self.dict_oxi
      #      dict_oxi[dop_el] = dop_dict['oxi_state']
      #      stp.add_oxidation_state_by_element(dict_oxi)
      #      E_ev = EwaldElectrostaticModel().get_energy(stp)
      #      ew.append(E_ev)


         #   st = Structure().update_from_pymatgen(stp)
         #   st_e.append((st, E_ev))      
      #   st_e = sorted(st_e, key=lambda list: list[1])
      #   st = []
      #   ew = []
      #   rr = []
      #   st_numb = min(st_numb, len(st_e))
      #   for i in range(st_numb):
      #    st.append(st_e[i][0])
      #    ew.append(st_e[i][1])
      #    rr
      #   ew = [e - min(ew) for e in ew]
      #   Nu = min(st_numb, len(st_d))

        
      #   return (st_min_e, st_max_r), (st_min_ed, st_max_rd)

    



def calc_conc(st, E, matrix, ALK = 'Li',  T = 1000, dif=1):
  """
  INPUT:
     st - siman structure
     E - decomposition energy(eV)
     matrix - matrix element
     dodant - dopant element
     dif - diference in oxi
  OUTPUT:
     conc (float) - concentration 
  """
  NA = len(st.get_elements_by_el_name(ALK)) + dif
  M = len(st.get_elements_by_el_name(matrix)) + 1

  
  k = 8.6*10**(-5)
  conc = np.sqrt(NA/M) * np.exp(-E/(2*k*T))
  
  return( conc )