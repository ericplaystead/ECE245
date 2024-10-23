import math
import sys

VSUBT = 0.026 #Thermal Voltage. 0.026 V
BC = float(10800000000000000000000000000000)
K=float(.0000862) # Boltzman constant
TEMP=float(600)
EG=float(1.12) # bandgap of Si, 1.12 eV
LILNSUBLILI = 1e10
LILQ = 1.6e-19
KSUBS=11.8  #dialectric current of Si. wiki says 11.68?
EPSILON_SUB_0=8.85e-14 #F/cm
ROOMTEMP_kT=0.026 #in eV

def muchgreater(a, b):
    #Return True if a is 20x greater than b
    if math.prod([a, 1/b]) > 20.0: return True
    else: return False

def iccofsi(temp=.01, element="si"):
    
    #print(type(temp))
    prod = math.sqrt(math.prod([BC,math.pow(temp, 3),math.exp(int(-1)*(EG/math.prod([K,temp])))]))
    print('Intrinsic carrier concentration of Si at {} K is {:e}'.format(temp, prod))
    return [prod,"{:e}".format(prod)]

    #final = math.sqrt()

def donorConcentration(NsubD:float, lilnsublili=LILNSUBLILI):
    n=0.0
    p=0.0
    print(lilnsublili)
    print("{:e}".format(lilnsublili))
    if muchgreater(NsubD, lilnsublili): n=NsubD
    else:
        print("NsubD: {} is not much greater than ni:{:e}".format(NsubD, lilnsublili))
        n = (NsubD/2) + math.pow((math.pow(NsubD,2)/4 + math.pow(lilnsublili, 2)), 0.5) 
    #print(math.pow(lilnsublili,2))
    p = (math.pow(lilnsublili,2))/n
    print("p: {:e}    n: {:e}".format(p,n))
    return (p,n)

def acceptorConcentration(NsubA:float, lilnsublili=LILNSUBLILI):
    n=0.0
    p=0.0
    if muchgreater(NsubA, lilnsublili): p=NsubA
    else:
        print("NsubA: {} is not much greater than ni:{:e}".format(NsubA, lilnsublili))
        p = (NsubA/2) + math.pow((math.pow(NsubA,2)/4 + math.pow(lilnsublili, 2)), 0.5)     
    n = (math.pow(lilnsublili,2))/p
    print("p: {:e}    n: {:e}".format(p,n))
    return (p,n)

def case3(NsubD: float, NsubA:float, lilnsublili=LILNSUBLILI):
    #case checking function for compensated semiconductors doped with donor and acceptor atoms
    p=0.0
    n=0.0

    if math.isclose(NsubD, lilnsublili, rel_tol=.0500) and math.isclose(NsubA, lilnsublili, rel_tol=.0500):
        print("NsubD and NsubA are too close to Nsubi, not n or p")
    elif muchgreater(NsubA - NsubD, lilnsublili):
        print("case3a:     p-type")
        if muchgreater(NsubA, NsubD):
            p=NsubA
            n=math.pow(lilnsublili, 2) / p
        else:
            p= NsubA-NsubD
    elif muchgreater(NsubD - NsubA, lilnsublili):
        print("case3b:     n-type")
        if muchgreater(NsubD, NsubA):
            n=NsubD
            p=math.pow(lilnsublili, 2) / n
        else: 
            p=NsubD-NsubA
    elif not muchgreater(math.fabs(NsubD - NsubA), lilnsublili):
        if NsubD > NsubA and NsubD > lilnsublili:
            print("case3c-1:     n-type")
            n = ((NsubD-NsubA) + math.sqrt( math.pow(NsubD-NsubA, 2)-(4 * math.pow(lilnsublili,2))))/2
            p = (math.pow(lilnsublili, 2)) / n
        elif NsubA > NsubD and NsubA > lilnsublili:
            print("case3c-2:     p-type")
            p = ((NsubA-NsubD) + math.sqrt( math.pow(NsubA-NsubD, 2)-(4 * math.pow(lilnsublili,2))))/2
            n = (math.pow(lilnsublili, 2)) / p

    print("p: {:e}    n: {:e}".format(p,n))
    return(p,n)
    
def electron_mobility(NsubT: float):
    # Only valid for Si at room temp
    mu_sub_n = 52.2+(1365/(1+math.pow((NsubT/(9.68*1e16)), 0.68)))    
    return mu_sub_n

def hole_mobility(NsubT: float):
    # Only valid for Si at room temp
    mu_sub_p = 44.9+(426/(1+math.pow((NsubT/(2.23*1e16)), 0.76)))    
    return mu_sub_p

def depletion_width(NsubA: float, NsubD: float, Vbi:float, VsubA: float):
    # Find depletion width of sides n and p, and also return total
    
    #return error if VsubA is greater than Vbi
    if math.fabs(VsubA)>Vbi:
        print("ERROR: applied voltage cannot be greater than built in voltage")
        return (0,0,0)
    xp = math.pow((((2*KSUBS*EPSILON_SUB_0)/LILQ)*(NsubD/(NsubA*(NsubA+NsubD)))*(Vbi-VsubA)), .5)
    
    xn = math.pow((((2*KSUBS*EPSILON_SUB_0)/LILQ)*(NsubA/(NsubD*(NsubA+NsubD)))*(Vbi-VsubA)), .5)
    WD = xp+xn
    return (xp, xn, WD)
    
def builtin_voltage(NsubA: float,NsubD: float):
    vbi = VSUBT * math.log((NsubA*NsubD)/math.pow(LILNSUBLILI,2))
    print("Calculated Built-in Voltage (Vbi): {:0.3f} V".format(vbi))
    return vbi

if __name__ == '__main__':
    
    NsubD=0.0
    NsubA=0.0
    NsubT=0.0
    temp=0.0

    print("""
 [1] - Intrinsic concentration calc
 [2] - Doped semiconductor calc
 [3] - Compensated semiconductor calc
 [4] - Electron mobility  (mu sub n)
 [5] - Hole mobility      (mu sub p)
 [6] - Built-in Potential V sub bi
 [7] - Depletion Region Width
        """)
    foo = input()
    if foo == "1":
        print("temp(K):     blank for room temp")
        t = input()
        if t:print(iccofsi(float(t)))
        else: print("ni at room temperature is 1e10 K")
    elif foo == "2":
        print("""
 [1] - Donor
 [2] - Acceptor""")
        atom = input()
        print("temp(K):    blank for room temp")
        temp = input()
        if temp: print("hitting on temp")
        if atom == "1":
            print("NsubD:")
            NsubD = float(input())          
            if temp: print(donorConcentration(NsubD, iccofsi(float(temp))[0]))
            else: print(donorConcentration(NsubD))

        elif atom == "2":
            print("NsubA:")
            NsubA = float(input())
            if temp: print(acceptorConcentration(NsubA, iccofsi(float(temp))[0]))
            else: print(acceptorConcentration(NsubA))

    elif foo == "3":
        print("NsubD:")
        NsubD = float(input())
        print("NsubA:")
        NsubA = float(input())
        print("temp(K):    blank for room temp")
        temp = input()
        if temp: print(case3(NsubD, NsubA, iccofsi(float(temp))[0]))
        else: print(case3(NsubD, NsubA))
        
    elif foo == "4":
        print("NsubT:")
        NsubT = float(input())
        print(electron_mobility(NsubT))
    elif foo == "5":
        print("NsubT:")
        NsubT = float(input())
        print("hole mobility mu sub p: {:0.5f}  cm^2 / V-s".format(hole_mobility(NsubT)))
    elif foo == "6":
        print("NsubD (n-side):")
        NsubD = float(input())
        print("NsubA (p-side):")
        NsubA = float(input())
        #print("temp(K):    blank for room temp")
        #temp = input()
        print("Calculated Built-in Voltage (Vbi): {:0.3f} V ".format(builtin_voltage(NsubA,NsubD)))
    elif foo == "7":
        print("NsubD (n-side):")
        NsubD = float(input())
        print("NsubA (p-side):")
        NsubA = float(input())
        #print("temp(K):    blank for room temp")
        #temp = input()
        print("doing room temp")
        print("Applied voltage (VA):    blank for 0 V")
        VsubA = input()
        if not VsubA: VsubA = 0.0
        foo=depletion_width(NsubA, NsubD, builtin_voltage(NsubA,NsubD),float(VsubA))
        print("""p-side depletion width(xp): {:0.3E}  
n-side depletion width(xn): {:0.3E}   
total depletion width(WD): {:0.3E}""".format(foo[0],foo[1],foo[2]))

    """try:
        ans=iccofsi(float(sys.argv[1]))
    except IndexError:
        ans=iccofsi()"""
    
    # Some examples of the functions.
    #print(math.isclose(NsubD, NsubA, rel_tol=.05))
    #print(NsubD)
    #print(NsubA)
    #print(math.prod([NsubA, 1/NsubD]))
    #print(math.prod([NsubD, 1/NsubA]))
    #print(muchgreater(NsubA, NsubD))
    #print(case3(NsubD, NsubA))
    #print(case3(NsubA, NsubD))
