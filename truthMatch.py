import sys, math
from array import array

from ROOT import TFile, TGraph, TCanvas, TPad, TLegend, TVector3, TH1D, TH2D
import ROOT
bec_border = 1450

def truthMatch(truth_positions, sp_positions, osp_positions, endCap = False):

    not_found, sfound, ofound = 0, 0, 0
    lx_found, ly_found, lx_ofound, ly_ofound, lx_nfound, ly_nfound = [], [], [], [], [], []

    truthPairs = find_pairs(truth_positions,12)
    diffs = []
    # print (len(truth_positions))
    for gpos1, l_pos in truth_positions:

        gpos2 = next((tpos2 if gpos1 == tpos1 else tpos1 for tpos1, tpos2 in truthPairs 
                      if gpos1 == tpos1 or gpos1 == tpos2), TVector3())
        
        matched_sp, sp_type = find_closest_match(gpos1, sp_positions, osp_positions,endCap)

        if matched_sp is None:
            diff = (999,999,999,999)            
            diffs.append(diff)
            continue
    
        dr = matched_sp.Perp() - gpos1.Perp()
        # dr = (matched_sp - gpos1).Perp()        
        delta = matched_sp - gpos1
        delta2 = matched_sp - gpos2

        if (delta.Perp() > delta2.Perp()):
            continue
            # delta = delta2

        diff = (delta.X(), delta.Y(), delta.Z(), dr)
        if endCap:
            # Assuming matched_sp and gpos1 have methods to return phi (azimuthal angle)
            dphi = matched_sp.Phi() - gpos1.Phi()
            # Normalize dphi to be in the range [-pi, pi] if necessary
            while dphi > math.pi:
                dphi -= 2*math.pi
            while dphi < -math.pi:
                dphi += 2*math.pi
            # print (diff)

            diff = (delta.X(), delta.Y(), dphi, dr)
        else:
            diff = (delta.X(), delta.Y(), delta.Z(), dr)
        
        diffs.append(diff)
        if (not endCap):
            if( abs(delta.Z()) > 12.0 and abs(delta.Z()) < 12.2):
                ofound += 1
                lx_ofound.append(l_pos[0])
                ly_ofound.append(l_pos[1])

            if( abs(delta.Z()) > 24.0 and abs(delta.Z()) < 24.4):
                ofound += 1
                lx_ofound.append(l_pos[0])
                ly_ofound.append(l_pos[1])

        if dr < 1:
            if sp_type == 'nominal':
                sfound += 1
                lx_found.append(l_pos[0])
                ly_found.append(l_pos[1])
            else:
                if(endCap):
                    ofound += 1
                    lx_ofound.append(l_pos[0])
                    ly_ofound.append(l_pos[1])
        else:
            not_found += 1
            lx_nfound.append(l_pos[0])
            ly_nfound.append(l_pos[1])

    return (
        sfound, ofound, not_found,
        lx_found, ly_found, lx_ofound, ly_ofound, 
        lx_nfound, ly_nfound,diffs )

def find_closest_match(position, sp_positions, osp_positions, endCap = False):
    min_diff = float('inf')
    matched_sp = None
    sp_type = ''

    for sp_pos in sp_positions + osp_positions:
        diff = (position - sp_pos).Perp()
        if (endCap):
            diff = abs((position - sp_pos).Z())
        if diff < min_diff:
            min_diff = diff
            matched_sp = sp_pos
            sp_type = 'nominal' if sp_pos in sp_positions else 'overlap'

    return matched_sp, sp_type


def compute_r(pos):
    return (pos.X()**2 + pos.Y()**2)**0.5

def split_positions(all_positions):
    if all_positions and isinstance(all_positions[0], tuple) and len(all_positions[0]) == 2:
        barrel = [(pos, lpos) for pos, lpos in all_positions if abs(pos.Z()) < bec_border]
        endcap = [(pos, lpos) for pos, lpos in all_positions if abs(pos.Z()) >= bec_border]
    else:
        barrel = [pos for pos in all_positions if abs(pos.Z()) < bec_border]
        endcap = [pos for pos in all_positions if abs(pos.Z()) >= bec_border]
    return barrel, endcap

def extractGlobalPos(positions):
    if positions and isinstance(positions[0], tuple) and len(positions[0]) == 2:
        x_values = [pos.X() for pos, _ in positions]
        y_values = [pos.Y() for pos, _ in positions]
        z_values = [pos.Z() for pos, _ in positions]
    else:
        x_values = [pos.X() for pos in positions]
        y_values = [pos.Y() for pos in positions]
        z_values = [pos.Z() for pos in positions]

    r_values = [math.sqrt(x*x + y*y) for x, y in zip(x_values, y_values)]
    return x_values, y_values, z_values, r_values

def extractLocalPos(positions):
    if positions and isinstance(positions[0], tuple) and len(positions[0]) == 2:
        x_values = [pos[0] for _,pos in positions]
        y_values = [pos[1] for _,pos,in positions]
        return x_values, y_values
def find_pairs(truthPositions, max_distance):
    distances = []
    for i in range(len(truthPositions)):
        for j in range(i + 1, len(truthPositions)):
            distance = calculate_distance(truthPositions[i][0], truthPositions[j][0])
            if distance < max_distance:
                distances.append((distance, i, j))
    
    distances.sort(key=lambda x: x[0])
    
    pairs = []
    used_indices = set()
    
    for _, i, j in distances:
        if i not in used_indices and j not in used_indices:
            pairs.append((truthPositions[i][0], truthPositions[j][0]))
            used_indices.add(i)
            used_indices.add(j)
    
    return pairs
def calculate_distance(pos1, pos2):
    diff_x = pos1[0] - pos2[0]
    diff_y = pos1[1] - pos2[1]
    distance = (diff_x**2 + diff_y**2 )**0.5
    return distance

