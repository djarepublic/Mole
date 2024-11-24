from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.analysis.energy_models import EwaldElectrostaticModel
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from siman.core.structure import Structure

def ewald(stp):
    Ew = EwaldElectrostaticModel()
    E = Ew.get_energy(stp)
    return(E)

# эта функция дает список неэквивалентных позиций определенног атома
def uniq_list(stp, atom):
    #stp = st.convert2pymatgen()
    list_uniq_all = SpacegroupAnalyzer(stp).get_symmetry_dataset()['equivalent_atoms']
    indices = [i for i, specie in enumerate(stp.species) if specie.symbol == atom]
    uniq_list = []
    for i in range(len(indices)):
        if indices[i] in list_uniq_all:
            uniq_list.append(indices[i])
    return (uniq_list)

def vac(stp, atom_v):
    #this list contain st and ew_energy
    #stp = st.convert2pymatgen()
    list_vac = []
    list_at_vac = uniq_list(stp, atom_v)
    for pos in list_at_vac:
        stp_c = stp.copy()
        st_new = stp_c.remove_sites([pos, ])
        list_vac.append(st_new)
 #   list_vac = sorted(list_vac, key=lambda list: list[1])
    return(list_vac)

def dop(stp, atom_d, atom_r):
    #this list contain st and ew_energy
    list_at_rep = uniq_list(stp, atom_r)
    list_dop = []
    for pos in list_at_rep:
        stp_c = stp.copy()
        st_new = stp_c.replace(pos, atom_d)
        list_dop.append(st_new)
    #list_dop = sorted(list_dop, key=lambda list: list[1])
    return (list_dop)

def dop_vac(st, atom_v, atom_d, atom_r):
    stp = st.convert2pymatgen()
    list_dop_vac = []
    list_dop = dop(stp, atom_d, atom_r)
    for st_dop in list_dop:
        st_dop_vac = vac(st_dop, atom_v)
        list_dop_vac = list_dop_vac + st_dop_vac

    list_dv_e = []
    for st_dop_vac in list_dop_vac:
        st_dv = Structure().update_from_pymatgen(st_dop_vac)
        st_dop_vac.add_oxidation_state_by_guess()
        E = ewald(st_dop_vac)
        list_dv_e.append((st_dv, E))
    list_dv_e = sorted(list_dv_e, key=lambda list: list[1])
    return(list_dv_e)


def dop2(st, atom_d1, atom_d2, atom_r1, atom_r2):
    stp = st.convert2pymatgen()
    list_dop1 = dop(stp, atom_d1, atom_r1)
    list_d1d2 = []
    for st_d1 in list_dop1:
        list_d1d2 = list_d1d2 + dop(st_d1, atom_d2, atom_r2)

    list_d1d2_e = []
    for stp_d1d2 in list_d1d2:
        st_d1d2 = Structure().update_from_pymatgen(stp_d1d2)
        stp_d1d2.add_oxidation_state_by_guess()
        E = ewald(stp_d1d2)
        list_d1d2_e.append((st_d1d2, E))
    list_d1d2_e = sorted(list_d1d2_e, key=lambda list: list[1])
    return(list_d1d2_e)


    
    



    


        
    