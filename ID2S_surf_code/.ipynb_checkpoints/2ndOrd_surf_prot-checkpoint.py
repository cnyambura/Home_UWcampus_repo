# Second order protein surf characterization

import MDAnalysis as mda
import statsmodels as stats
import math
import numpy as np
import pandas as pd


def solv_vec_protsurf(protein_surf, allprot_atoms, near_neighb, rsear_1, rsear_2, rd_space):
    
    """Calculate solvent vector for surface residue in AtomGroup protein_surf, as described in
    Jones, S. and J. M. Thornton (1997). "Analysis of protein-protein interaction sites using 
    surface patches." J Mol Biol 272(1): 121-132. 
    
    Here, the solvent vector is calculated by taking the a single central C alpha atom of a surface
    residue and its nearest atom neighbors and subsquently calulating the center of mass. The inverse 
    of the vector from the C alpha atom of the central surface residue to the calc. center of mass was 
    then calculated. 
    
    Inputs consisted of: 
    
    - protein_surf: Atom group from MD Analysis containing only the protein surface atoms 
    - allprot_atoms: Atom group from MD Analysis containing all protein atoms  
    - near_neighb: Desired umber of nearest neighbors for each surface residue C alpha atom
    - rsear_1: Lower radius bound (in Angstroms) for MDAnalysis Neighbor Search to find # of nearest neighbors
    - rsear_2: Upper radius bound (in Angstroms) for MDAnalysis Neighbor Search to find # of nearest neighbors
    - rd_space: Spacing value for the radii array used in the neighbor search 
    
    Output consists of: 
    
    - new_slv: dictionary containing the solvent vector(value) for each surface residue(key)
    - nn_res: Actual number of nearest neighbors used for the vector calc., for each residue 
    
    """
    
    # Import MDAnalysis 
    import MDAnalysis as mda
    
    #Import numpy 
    import numpy as np
    
    #Initialize dictionary 
    new_slv = {}
    
    # Get number of surface residues
    ln_surfres = len(list(protein_surf.residues))
    
    # list to store nearest neighbors for each residue used for the nearest neighbor search 
    nn_res = []
    
    # For each surface residue, we need to calculate the solvetn vector
    for i in range(ln_surfres):
        
        # select current residue 
        resAtms = protein_surf.residues[i].atoms
        
        # Select C alpha atom
        CA_res = resAtms.select_atoms('name CA')
        
        # Use Atom Neighbor search to get its 10 nearest neighbors
        
        # Initialize AtomNeighbor Search class
        nn_psurf = mda.lib.NeighborSearch.AtomNeighborSearch(allprot_atoms)
        
        # run for loop to find radius that will give me 10 nearest neighbors
        
        # Initialize numpy array of even spaced radius values in Angstroms
        rs_arr = np.arange(rsear_1, rsear_2, rd_space)
        rad_indx = 0
        
        # Looping through radius values, save index that points to radius value that give user input NN
        for j in range(len(rs_arr)):
            
            chk_ca = list(nn_psurf.search(atoms=CA_res, radius=rs_arr[j], level='A'))
            
            # I need to add elif 11 neighbors because based on the cutoff, some surf res will have 
            # 9,11, or 12 neighbors but not 10 neighbors, Can't think of a better fix
            if near_neighb <= len(chk_ca) <= near_neighb+2:
                
                rad_indx = j
                #print(len(chk_ca))
                nn_res.append(len(chk_ca))
                
                break 
        
        # Error check: Making sure that a radial value was stored that corresponds to desired nearest neighbors
        if rad_indx == 0:    
            print(rad_indx)
            print(list(resAtms))
            print(i)
        
        # Atom group containing NN around the C alpha atom 
        neighb_fd = nn_psurf.search(atoms=CA_res, radius=rs_arr[rad_indx], level='A')
        
        # Solvent Vector is the inverse of COM - CA vector or CA - COM vector (or (COM - CA)*-1)
        slv = CA_res.positions - neighb_fd.center_of_mass()
        
        slv_norm = slv/np.linalg.norm(slv)
        
        #print(slv.shape)
        
        if slv.shape == (2,3):
            print(slv.shape)
            print(CA_res)
        
        
        if np.all(slv == 0) == True:
            print(neighb_fd)
            
        # Save solvent vector in dictionary 
        new_slv[str(protein_surf.residues[i])] = slv_norm
        
    return new_slv, nn_res