
from pymatgen.analysis.local_env import  CrystalNN
from siman.geo import supercell
import numpy as np 
from pymatgen.io.ase import AseAtomsAdaptor
from .mat_data import  eln, exceptions, table_ir, rn, ionic_radius
from pymatgen.io.vasp import Chgcar
from mp_api.client import MPRester
import os

from pymatgen.analysis.defects.generators import ChargeInterstitialGenerator

from siman.core.structure import Structure
pymatgen2siman = Structure().update_from_pymatgen


class Dopant_dv:
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
        
        dict_oxi = {}
        for i, el in enumerate(st.get_elements()):
            dict_oxi[el] = oxi_list[i]
        self.dict_oxi = dict_oxi

        self.st = st
        self.ion_types = {}

        for at in at_list:
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
       
       

    def search_dop(self, prec_r, prec_eln, dif=1, reduce = 2):   
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
                                (prop['r_ionic']>r_min and prop['r_ionic']<r_max) and
                                (abs(m_eln - d_eln) < prec_eln) and (d_el not in exceptions) and
                                (d_el not in self.ion_types) 
                               ):
                              
                              if( 
                                 ((dif==1) and (self.ion_types[ion]['oxi_state'] +1 == d_oxi )) or 
                                  ((dif==2) and (self.ion_types[ion]['oxi_state'] +2 == d_oxi )) or
                                  ((self.ion_types[ion]['oxi_state'] +1 == d_oxi or self.ion_types[ion]['oxi_state'] +2 == d_oxi) and (dif==None))
                                 ):
                               
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

         if reduce == 2:
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
            
    def make_dop(self, dop_el, dop_dict):
        """
        INPUT:
           dop_el (str) - dopant element
           dop_dict (dict) - dopant_dict from search_dop    
        OUTPUT:
           std - structure with substituted dopant atom

        """
        
        
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
               

        in_n = c_st.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
        for j in range(len(in_n)):
           if in_n[j]==dop_dict['dop_position']:
              new_pos = j
              break
           
        std = c_st.replace_atoms([new_pos], dop_el) 
        return std

       

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
        r = []
        xred = []

        cst = self.make_dop(dop_el, dop_dict)   
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


    def mine_st(self, dop_el, dop_dict, calc = None):
            """
             INPUT:
                dop_el (str) - dopant element
                dop_dict (dict) - dopant_dict from search_dop
                draw - add Pu atom in vacancy
    
             OUTPUT:
                st_min_e, st_min_ed - structure (dopant+vacancy) with minimum energy, (st_min_ed for visualization)
    
            """
            
           
            r_sorted, st_sorted, xred_sorted = self.make_st(dop_el, dop_dict)
            st_short_mine = st_sorted[:4]
            st_dp_mine = [st.convert2pymatgen() for st in st_short_mine]
            st_da_mine = [AseAtomsAdaptor.get_atoms(s) for s in st_dp_mine]

            # sevennet_0_cal = SevenNetCalculator("7net-0") 

            E_UP_mine = []
            for s in st_da_mine:
               s.calc = calc
               E_UP_mine.append(s.get_potential_energy())

            ind_min_e = E_UP_mine.index(min(E_UP_mine))
            st_min_e = st_sorted[ind_min_e]
            st_min_ed = st_min_e.add_atom(xred_sorted[ind_min_e ], 'Pu')            
            return st_min_e, st_min_ed
    


    def maxr_st(self, dop_el, dop_dict, draw = False):
       """
        INPUT:
           dop_el (str) - dopant element
           dop_dict (dict) - dopant_dict from search_dop
           draw - add Pu atom in vacancy
 
        OUTPUT:
           st_max_r, st_max_rd - structure (dopant+vacancy) with maximum dopant-vacancy distance, (st_max_rd for visualization)
 
       """
       
       r_sorted, st_sorted, xred_sorted = self.make_st(dop_el, dop_dict)
       ind_max = len(st_sorted) - 1
       st_max_r = st_sorted[ind_max ]
       st_max_rd = st_max_r.add_atom(xred_sorted[ind_max ], 'Pu')
       return st_max_r, st_max_rd
           
    


class Defects_di:
    def __init__(self, alk, id=None, chgcar_file = None):
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
        try:
            # Access the environment variable named 'API_KEY'
            MP_API_KEY = os.environ['MP_API_KEY']
        except KeyError:
            print("Error: MP_API_KEY environment variable not set.")

        if id:
            with MPRester(MP_API_KEY) as mpr:
              chgcar = mpr.get_charge_density_from_material_id(id)
         # chgcar.write_file("CHGCAR_"+id)
        if chgcar_file:
           chgcar = Chgcar.from_file(chgcar)
           
        stp =  chgcar.structure
        st = pymatgen2siman(stp)

        self.st = st.copy() 
        self.alk = alk
      #   st = st.get_primitive_cell()
        at_list = st.get_unique_type_els()
        oxi_list = st.get_oxi_states('guess')
      #   stp = st.convert2pymatgen()
        
        dict_oxi = {}
        for i, el in enumerate(st.get_elements()):
            dict_oxi[el] = oxi_list[i]
        self.dict_oxi = dict_oxi

        self.ion_types = {}

        for at in at_list:
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

        cig = ChargeInterstitialGenerator(max_insertions = 1)
        ins = cig.generate(chgcar, insert_species=["Pu"])
        defect = next(ins)
        defect.site.coords
         #     new_structure = chgcar.structure.copy()
        self.sti = self.st.add_atom(defect.site.frac_coords, 'Pu')

            

    def look_match(self):
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
       
   #  def make_inst(self, id, chgcar):

         #  new_structure.to('POSCAR1')
          
          

    def search_dop(self, prec_r, prec_eln, dif=1, reduce = 1):   
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
                                (prop['r_ionic']>r_min and prop['r_ionic']<r_max) and
                                (abs(m_eln - d_eln) < prec_eln) and (d_el not in exceptions) and
                                (d_el not in self.ion_types) 
                               ):
                              
                              if( 
                                 ((dif==1) and (self.ion_types[ion]['oxi_state']  == d_oxi + 1 )) or 
                                  ((dif==2) and (self.ion_types[ion]['oxi_state'] == d_oxi  + 2)) or
                                  ((self.ion_types[ion]['oxi_state']  == d_oxi+1 or self.ion_types[ion]['oxi_state']  == d_oxi + 2) and (dif==None))
                                 ):
                               
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
         
       if reduce == 2:
         dop2 = []
         try:
            rr = self.look_match()
            if rr:
               if rr[0]< 0.05:
                  target_key = 'oxi_state'
                  compare = [rr[1][0], rr[1][1]]
                  #print(compare)
                  atom = max(compare, key=lambda k: self.ion_types[k][target_key])
                  print(atom)
                  for item in dop:
                     if item['matrix_el'] != atom:
                        dop2.append(item)
                  dop = dop2
         except:
            pass
         
       
       return(dop)
            
    def make_dop(self, dop_el, dop_dict):
        """
        INPUT:
           dop_el (str) - dopant element
           dop_dict (dict) - dopant_dict from search_dop    
        OUTPUT:
           st_minr - structure with substituted dopant atom and alk interstitial with minimum distance
           st_maxr -  maximum distance
           st_minrd -  same as st_minr inst Li replaced by Pu for visualization
           st_maxrd - same as st_maxr inst Li replaced by Pu for visualization

        """
        
        
        if max(self.sti.rprimd_len()) < 10:
            l = [10, 10, 10]
            c_st = supercell(self.sti, l, test_natom = 0, skip_numbers = [self.sti.natom - 1])
        else:
            l = []
            for i in range(3):
               if self.sti.rprimd_len()[i]<10:
                  l.append(10)
               else:
                  l.append(self.sti.rprimd_len()[i])
            c_st = supercell(self.sti, l, test_natom = 0, skip_numbers = [self.sti.natom - 1])
        in_list = c_st.get_el('Pu')
               

        in_n = c_st.init_numbers        # эотот блок востанавливает исходные номера новой структуры по старой
        dop_pos = []
        for j in range(len(in_n)):
           if in_n[j]==dop_dict['dop_position']:
              dop_pos.append(j)
         
        r_min = 100
        r_max = 0
        dpmin = None
        dpmax = None
        for dp in dop_pos:
           r = c_st.distance(in_list[0], dp)
           if r>r_max:
              r_max = r
              dpmax = dp
           if r<r_min:
              r_min = r
              dpmin = dp
 
           
        st_minrd = c_st.replace_atoms([dpmin], dop_el)
        st_maxrd = c_st.replace_atoms([dpmax], dop_el)
        
        inst_pos = st_minrd.get_el('Pu')[0]
        st_minr = st_minrd.replace_atoms([inst_pos], self.alk )

        inst_pos = st_maxrd.get_el('Pu')[0]
        st_maxr = st_maxrd.replace_atoms([inst_pos], self.alk )
        
      #   st_minr = st_minr.remove_atoms(in_list[1:])
      #   st_maxr = st_maxr.remove_atoms(in_list[1:])
         
        return st_minr, st_maxr, st_minrd, st_maxrd



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

def Cl(ch, Eb, T=300, g = 6):
    """
    INPUT:
       ch - hight temperature vacancy concentration (c_vac = N_matrix/N_dop * c_dop)
       Eb - binding energy(eV)
    OUTPUT:
       conc (float) - free RT vacancy concentration 
    """
    k = 8.617e-05
    al = np.exp(-Eb/(k*T))/g
    C = (np.sqrt(al**2 +4*al*ch) - al) / (2 * (1 - al))
    return C
