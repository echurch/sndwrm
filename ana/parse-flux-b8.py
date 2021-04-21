
#curl -G http://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt > ../macros/b8.dat
#curl -G http://www.sns.ias.edu/~jnb/SNdata/Export/hepspectrum/hepspectrum.txt > ../macros/hep.dat


from ROOT import TFile, TH1F
import pdb

hname = "b8"
hout = TH1F(hname,"b8-jb",50,0.0,16.0)
fn = "../macros/b8.dat"
data = []
energy = []
dots = False

with open(fn,'r') as f:
    lines = f.readlines()
    for line in lines:
        if "....." not in line and not dots:
            continue

        dots = True

        if line == "\n":
            continue


        if "...." not in line and len(line.split()) == 4:
            energy.append(float(line.split()[0]))
            data.append(float(line.split()[1]))

print("Filling histogram, saving to file. Also dumping for cnp'ing into MARLEY's json file, since somehow it's not working w ROOT.")
for ii in range(len(energy)):
    hout.Fill(energy[ii],data[ii])
print("ernergy, data lengths" + str(len(energy)) + ", " +str(len(data)))
print(str(energy[::8]),sep=",")
print(str(data[::8]),sep=",")


fon = "../macros/b8_bahcall.root"
fout = TFile(fon,"recreate")
hout.Write()
fout.Close()

# Read it back and check.

fin = TFile(fon)

# form fin.Get("b8").GetNbinsX()
hreadfunc = "fin" + "." + "Get(" + "'" + hname + "'" + ")" + "." 

exit()

## fin.b8.Draw() confirms all's well.
pdb.set_trace()  
for bin in range(eval(hreadfunc + "GetNbinsX()")):
    print ("bin/data: " + str(eval(hreadfunc + "GetBin(ii)")) + "/" + str(eval(hreadfunc + "GetBinContent(ii)")) )


