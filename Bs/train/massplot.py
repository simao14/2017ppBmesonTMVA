import ROOT as r


if 'can1' not in locals():
    can1 = r.TCanvas('can1', 'canvas', 800, 600)
r.gStyle.SetOptStat(0)
r.gStyle.SetHistMinimumZero()

fin = r.TFile('massplot.root')

bmass = 5.36682
pts = [7, 10, 15, 20, 30, 50]
# pts = [7, 10]
for i in range(len(pts) - 1):
    hname = f'massa_{pts[i]}_{pts[i+1]}'
    h = fin.Get(hname)
    h.SetLineWidth(2)
    h.SetXTitle('m_{B} [GeV]')
    h.Draw()
    latex = r.TLatex()
    latex.DrawLatexNDC(0.2, 0.8, f'{pts[i]} < ' + r'p_{T}' + f'< {pts[i + 1]}')
    xli = h.GetXaxis().FindBin(bmass - 0.3)
    xlf = h.GetXaxis().FindBin(bmass - 0.2)
    xri = h.GetXaxis().FindBin(bmass + 0.2)
    xrf = h.GetXaxis().FindBin(bmass + 0.3)
    nside = h.Integral(xli, xlf)
    nside += h.Integral(xri, xrf)
    print(nside)
    latex.DrawLatexNDC(0.6, 0.8, f'{nside:.0f} events ')
    latex.DrawLatexNDC(0.6, 0.74, 'in sideband')
    latex.DrawLatexNDC(0.6, 0.68, 'after preselection')
    can1.SaveAs(f'aftermass{pts[i]}_{pts[i+1]}.png')

for i in range(len(pts) - 1):
    hname = f'massb_{pts[i]}_{pts[i+1]}'
    h = fin.Get(hname)
    h.SetLineWidth(2)
    h.SetXTitle('m_{B} [GeV]')
    h.Draw()
    latex = r.TLatex()
    latex.DrawLatexNDC(0.2, 0.8, f'{pts[i]} < ' + r'p_{T}' + f'< {pts[i + 1]}')
    xli = h.GetXaxis().FindBin(bmass - 0.3)
    xlf = h.GetXaxis().FindBin(bmass - 0.2)
    xri = h.GetXaxis().FindBin(bmass + 0.2)
    xrf = h.GetXaxis().FindBin(bmass + 0.3)
    nside = h.Integral(xli, xlf)
    nside += h.Integral(xri, xrf)
    print(nside)
    latex.DrawLatexNDC(0.6, 0.8, f'{nside:.0f} events ')
    latex.DrawLatexNDC(0.6, 0.74, 'in sideband')
    latex.DrawLatexNDC(0.6, 0.68, 'before preselection')
    can1.SaveAs(f'beforemass{pts[i]}_{pts[i+1]}.png')

# can1.Update()
