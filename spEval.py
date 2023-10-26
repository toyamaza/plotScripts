#!/usr/bin/env python3

import sys, math
from array import array

from ROOT import TFile, TGraph, TCanvas, TPad, TLegend, TVector3, TH1D, TH2D
import ROOT
from truthMatch import *

bec_border = 1450

def TPGraph(n, x, y):
    xx = array('d', x)
    yy = array('d', y)
    return TGraph(n, xx, yy)


def getPositions(gl_x, gl_y, gl_z, loc_x=None, loc_y=None):
    positions = []
    for j in range(len(gl_x)):
        x, y, z = gl_x[j], gl_y[j], gl_z[j]
        gpos = TVector3(x, y, z)
        
        if loc_x is not None and loc_y is not None:
            lpos = (loc_x[j], loc_y[j])
            positions.append((gpos, lpos))
        else:
            positions.append(gpos)
            
    return positions

def create_histogram(data, nbins, xlow, xup, title):
    hist = ROOT.TH1D(title, title, nbins, xlow, xup)
    for value in data:
        hist.Fill(value)
    return hist

def draw_and_save_hist(hist, frame_range, x_title, y_title, filename):
    canv = ROOT.TCanvas("canv", "canv", 800, 600)
    frame = ROOT.gPad.DrawFrame(*frame_range)
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)
    color = ROOT.kRed if 'B' in filename else ROOT.kBlue
    hist.SetLineColor(color)
    hist.Draw("same")
    canv.SaveAs(filename)
    canv.Clear()

def extract_and_draw(diffs, name_prefix):

    dxs = [diff[0] for diff in diffs]
    dys = [diff[1] for diff in diffs]
    dzs = [diff[2] for diff in diffs]
    drs = [diff[3] for diff in diffs]

    xmax = 10
    xmax_big = 40
    # Create histograms
    hist_dx = create_histogram(dxs, 100, -xmax, xmax, "dx_" + name_prefix)
    hist_dy = create_histogram(dys, 100, -xmax, xmax, "dy_" + name_prefix)
    hist_dz = create_histogram(dzs, 1000, -xmax_big, xmax_big, "dz_" + name_prefix)
    hist_dr = create_histogram(drs, 100, -xmax, xmax, "dr_" + name_prefix)

    # Draw and save histograms
    draw_and_save_hist(hist_dx, (-xmax, 0, xmax, hist_dx.GetMaximum()*1.4 ), "dx", "Entries", "fig/res" + name_prefix + "_x.pdf")
    draw_and_save_hist(hist_dy, (-xmax, 0, xmax, hist_dy.GetMaximum()*1.4 ), "dy", "Entries", "fig/res" + name_prefix + "_y.pdf")
    draw_and_save_hist(hist_dz, (-xmax_big, 0, xmax_big, hist_dz.GetMaximum()*1.4 ), "dz", "Entries", "fig/res" + name_prefix + "_z.pdf")
    draw_and_save_hist(hist_dz, (-xmax_big, 0, xmax_big, hist_dz.GetMaximum()*1.4 ), "dz", "Entries", "fig/res" + name_prefix + "_z.root")    
    draw_and_save_hist(hist_dr, (-xmax, 0, xmax, hist_dr.GetMaximum()*1.4 ), "dr", "Entries", "fig/res" + name_prefix + "_r.pdf")




def create_tgraphs(x_values, y_values, z_values, r_values):
    zr_graph = TPGraph(len(z_values), z_values, r_values)
    xy_graph = TPGraph(len(x_values), x_values, y_values)
    return zr_graph, xy_graph

def create_tgraph(x_values, y_values):
    xy_graph = TPGraph(len(x_values), x_values, y_values)
    return xy_graph


    


def print_stripCluster_data(ntuple_file):

    ROOT.gROOT.LoadMacro("atlasstyle/AtlasStyle.C")
    ROOT.gROOT.LoadMacro("atlasstyle/AtlasLabels.C")
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 10001;")

    from ROOT import SetAtlasStyle
    SetAtlasStyle()
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetMarkerSize(0.2)

    f = TFile.Open(ntuple_file)
    tree = f.Get("TRKTree")
    n_entries = tree.GetEntries()

    grRecoClusterx, grRecoClustery, grRecoSPx, grRecoSPy = [], [], [], []
    grPairDx, grPairDy, grPairDz, grPairDr = [], [], [], []
    allTruthPositions, allClusterPositions, allSPPositions, allOSPPositions = [], [], [], []

    lx_foundE, ly_foundE, lx_ofoundE, ly_ofoundE, lx_nfoundE, ly_nfoundE = [], [], [], [], [], []
    foundE, oFoundE, notFoundE = 0, 0, 0

    lx_foundB, ly_foundB, lx_ofoundB, ly_ofoundB, lx_nfoundB, ly_nfoundB = [], [], [], [], [], []
    foundB, oFoundB, notFoundB = 0, 0, 0
    diffsE, diffsB = [], []
    for i in range(n_entries):
        tree.GetEntry(i)  # Load the i-th entry into memory

        # truth hit positions
        truthPositions = getPositions(
            tree.stripCluster_truth_globX,
            tree.stripCluster_truth_globY,
            tree.stripCluster_truth_globZ,
            tree.stripCluster_truth_locX,
            tree.stripCluster_truth_locY)

        # reco clusters
        clusterPositions =  getPositions(
            tree.stripCluster_reco_globX, 
            tree.stripCluster_reco_globY, 
            tree.stripCluster_reco_globZ,
            tree.stripCluster_reco_locX, 
            tree.stripCluster_reco_locY)

        # nominal SP
        spPositions = getPositions(
            tree.stripSpacePoint_reco_globX, 
            tree.stripSpacePoint_reco_globY, 
            tree.stripSpacePoint_reco_globZ)

        # overlap SP
        ospPositions = getPositions(
            tree.ostripSpacePoint_reco_globX, 
            tree.ostripSpacePoint_reco_globY, 
            tree.ostripSpacePoint_reco_globZ)
            
        truthPositionsB = [(pos, lpos) for pos, lpos in truthPositions if abs(pos.Z()) < bec_border]
        truthPositionsE = [(pos, lpos) for pos, lpos in truthPositions if abs(pos.Z()) >= bec_border]
        # print (len(truthPositions),len(truthPositionsB),len(truthPositionsE))
        resultsE = truthMatch(truthPositionsE, spPositions, ospPositions, True)
        foundE+=resultsE[0]
        oFoundE+=resultsE[1]
        notFoundE+=resultsE[2]
        lx_foundE.extend(resultsE[3])
        ly_foundE.extend(resultsE[4])
        lx_ofoundE.extend(resultsE[5])
        ly_ofoundE.extend(resultsE[6])

        lx_nfoundE.extend(resultsE[7])
        ly_nfoundE.extend(resultsE[8])


        # diffsE.append(resultsE[9])
        diffsE.extend(resultsE[9])        
        
        resultsB = truthMatch(truthPositionsB, spPositions, ospPositions)
        foundB+=resultsB[0]
        oFoundB+=resultsB[1]
        notFoundB+=resultsB[2]

        lx_foundB.extend(resultsB[3])
        ly_foundB.extend(resultsB[4])

        lx_ofoundB.extend(resultsB[5])
        ly_ofoundB.extend(resultsB[6])

        lx_nfoundB.extend(resultsB[7])
        ly_nfoundB.extend(resultsB[8])
        
        diffsB.extend(resultsB[9])

        allTruthPositions.extend(truthPositions)
        allClusterPositions.extend(clusterPositions)
        allSPPositions.extend(spPositions)
        allOSPPositions.extend(ospPositions)

        
    extract_and_draw(diffsE, "E")
    extract_and_draw(diffsB, "B")        

    print ("B eff = ",(foundB + oFoundB)/float(foundB+oFoundB + notFoundB), " " ,foundB, "( " ,foundB, ", ",oFoundB," )",notFoundB)    
    print ("E eff = ",(foundE + oFoundE)/float(foundE+oFoundE + notFoundE), " " ,foundE, "( " ,foundE, ", ",oFoundE," )",notFoundE)

    allTruthPositionsB, allTruthPositionsE = split_positions(allTruthPositions)
    allClusterPositionsB, allClusterPositionsE = split_positions(allClusterPositions)
    allSPPositionsB, allSPPositionsE = split_positions(allSPPositions)
    allOSPPositionsB, allOSPPositionsE = split_positions(allOSPPositions)

    arrTruthB_x, arrTruthB_y, arrTruthB_z, arrTruthB_r = extractGlobalPos(allTruthPositionsB)
    arrTruthE_x, arrTruthE_y, arrTruthE_z, arrTruthE_r = extractGlobalPos(allTruthPositionsE)

    arrTruthB_lx, arrTruthB_ly = extractLocalPos(allTruthPositionsB)
    arrTruthE_lx, arrTruthE_ly = extractLocalPos(allTruthPositionsE)
    
    arrClusterB_x, arrClusterB_y, arrClusterB_z, arrClusterB_r = extractGlobalPos(allClusterPositionsB)
    arrClusterE_x, arrClusterE_y, arrClusterE_z, arrClusterE_r = extractGlobalPos(allClusterPositionsE)

    arrSPB_x, arrSPB_y, arrSPB_z, arrSPB_r = extractGlobalPos(allSPPositionsB)
    arrSPE_x, arrSPE_y, arrSPE_z, arrSPE_r = extractGlobalPos(allSPPositionsE)

    arrOSPB_x, arrOSPB_y, arrOSPB_z, arrOSPB_r = extractGlobalPos(allOSPPositionsB)
    arrOSPE_x, arrOSPE_y, arrOSPE_z, arrOSPE_r = extractGlobalPos(allOSPPositionsE)

    tgTruthB_zr, tgTruthB_xy = create_tgraphs(arrTruthB_x, arrTruthB_y, arrTruthB_z, arrTruthB_r)
    tgTruthE_zr, tgTruthE_xy = create_tgraphs(arrTruthE_x, arrTruthE_y, arrTruthE_z, arrTruthE_r)

    tgTruthB_local_xy = create_tgraph(arrTruthB_lx, arrTruthB_ly)
    tgTruthE_local_xy = create_tgraph(arrTruthE_lx, arrTruthE_ly)
    # print (len(lx_foundB))
    # print (len(ly_foundB))
    # for i in range (len(lx_foundB)):
    #     print (lx_foundB[i], ly_foundB[i])
    tgTruthB_localFound = create_tgraph(lx_foundB,ly_foundB)
    tgTruthB_localOFound = create_tgraph(lx_ofoundB,ly_ofoundB)
    tgTruthB_localNFound = create_tgraph(lx_nfoundB,ly_nfoundB)        

    tgTruthE_localFound = create_tgraph(lx_foundE,ly_foundE)
    tgTruthE_localOFound = create_tgraph(lx_ofoundE,ly_ofoundE)
    tgTruthE_localNFound = create_tgraph(lx_nfoundE,ly_nfoundE)        

    tgClusterB_zr, tgClusterB_xy = create_tgraphs(arrClusterB_x, arrClusterB_y, arrClusterB_z, arrClusterB_r)
    tgClusterE_zr, tgClusterE_xy = create_tgraphs(arrClusterE_x, arrClusterE_y, arrClusterE_z, arrClusterE_r)

    tgSPB_zr, tgSPB_xy = create_tgraphs(arrSPB_x, arrSPB_y, arrSPB_z, arrSPB_r)
    tgSPE_zr, tgSPE_xy = create_tgraphs(arrSPE_x, arrSPE_y, arrSPE_z, arrSPE_r)

    tgOSPB_zr, tgOSPB_xy = create_tgraphs(arrOSPB_x, arrOSPB_y, arrOSPB_z, arrOSPB_r)
    tgOSPE_zr, tgOSPE_xy = create_tgraphs(arrOSPE_x, arrOSPE_y, arrOSPE_z, arrOSPE_r)

    zr_range = (-3000,0,3000,1500)
    xy_range = (-1200, -1200, 1200, 1200)


    
    canv = ROOT.TCanvas('c1','Graphs',800,600)


    def draw_and_save(frame_range, x_title, y_title, graph, filename):
        frame = ROOT.gPad.DrawFrame(*frame_range)
        frame.GetXaxis().SetTitle(x_title)
        frame.GetYaxis().SetTitle(y_title)
        color = ROOT.kRed if 'B' in filename else ROOT.kBlue
        graph.SetMarkerColor(color)        
        graph.Draw("p same")
        canv.SaveAs(filename)
        canv.Clear()
        
    def draw_and_save2(frame_range, x_title, y_title, graphB, graphE, filename):
        frame = ROOT.gPad.DrawFrame(*frame_range)
        frame.GetXaxis().SetTitle(x_title)
        frame.GetYaxis().SetTitle(y_title)

        graphB.SetMarkerColor(ROOT.kRed)  # Setting color for graph B
        graphB.Draw("p same")

        graphE.SetMarkerColor(ROOT.kBlue)  # Setting color for graph E
        graphE.Draw("p same")

        canv.SaveAs(filename)
        canv.Clear()

    def draw_and_save3(frame_range, x_title, y_title, graph1, graph2, graph3, filename):
        canv.Clear()
        canv.cd()
        frame = ROOT.gPad.DrawFrame(*frame_range)
        frame.GetXaxis().SetTitle(x_title)
        frame.GetYaxis().SetTitle(y_title)

        
        graph1.SetMarkerColor(ROOT.kRed)
        graph1.Draw("p same")


        graph2.SetMarkerColor(ROOT.kBlue)
        graph2.Draw("p same")

        graph3.SetMarkerColor(ROOT.kGreen)
        graph3.Draw("p same")

        canv.SaveAs(filename)



    graphs_zr = [
        (tgTruthB_zr,tgTruthE_zr, "fig/truth_zr.pdf"),
        # (tgClusterB_zr,tgClusterE_zr, "fig/clusterB_zr.pdf"),
        (tgSPB_zr,tgSPE_zr, "fig/SP_zr.pdf"),
        (tgOSPB_zr,tgOSPE_zr, "fig/OSP_zr.pdf"),        
    ]

    graphs_xy = [
        (tgTruthB_xy, "fig/truthB_xy.pdf"),
        (tgTruthE_xy, "fig/truthE_xy.pdf"),
        # (tgClusterE_xy, "fig/clusterE_xy.pdf"),
        (tgSPB_xy, "fig/SPB_xy.pdf"),
        (tgSPE_xy, "fig/SPE_xy.pdf"),
        (tgOSPB_xy, "fig/OSPB_xy.pdf"),
        (tgOSPE_xy, "fig/OSPE_xy.pdf")
    ]

    graphs_local = [
        (tgTruthB_local_xy, "fig/truthB_loc_xy.pdf"),
        (tgTruthE_local_xy, "fig/truthE_loc_xy.pdf"),
    ]
    
    # Draw and save for zr range
    for graphB, graphE, filename in graphs_zr:
        draw_and_save2(zr_range, "Z", "R", graphB, graphE, filename)

    # Draw and save for xy range
    for graph, filename in graphs_xy:
        draw_and_save(xy_range, "X", "Y", graph, filename)

    # for graph, filename in graphs_local:
    local_range = (-80, -80, 80, 80)        
    draw_and_save(local_range, "X", "Y", tgTruthB_local_xy, "fig/truthB_loc_xy.pdf")
    local_range = (-200, 300, 200, 1000)        
    draw_and_save(local_range, "X", "Y", tgTruthE_local_xy, "fig/truthE_loc_xy.pdf")
    local_range = (-80, -80, 80, 80)        
    draw_and_save3(local_range, "X", "Y", tgTruthB_localFound,tgTruthB_localOFound, tgTruthB_localNFound,"fig/truthB_localEff.pdf")
    local_range = (-200, 300, 200, 1000)            
    draw_and_save3(local_range, "X", "Y", tgTruthE_localFound,tgTruthE_localOFound, tgTruthE_localNFound,"fig/truthE_localEff.pdf")
    
    return


def main(argv):
    
    ntuple = 'SPNtuple_local2.root'
    # ntuple = 'SPNtuple.root'        
    print_stripCluster_data(ntuple)
    
if __name__ == '__main__':
    main(sys.argv)
