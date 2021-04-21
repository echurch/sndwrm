

from ROOT import TFile, TH1

#neutrons
#f0 = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_0.root")
f0 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_0.root")
h0ES = f0.H22
h0CC = f0.H23
#f1 = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_1.root")
f1 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_1.root")
h1ES = f1.H22
h1CC = f1.H23
#f2 = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_2.root")
f2 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_2.root")
h2ES = f2.H22
h2CC = f2.H23
#f3 = TFile("/Volumes/EC2TB/G4/mod4/neutron_coldcryoskin_3.root")
f3 = TFile("/Volumes/EC2TB/G4/mod4/marley_optphys_hep_3.root")
h3ES = f3.H22
h3CC = f3.H23


h0ES.Add(h1ES)
h0ES.Add(h2ES)
h0ES.Add(h3ES)

h0CC.Add(h1CC)
h0CC.Add(h2CC)
h0CC.Add(h3CC)

h0CC.Draw()
h0ES.Draw('s')
