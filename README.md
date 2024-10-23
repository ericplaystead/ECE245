# ECE245
A basic calculator and notes for BSU course ECE245 - Intro to Electronic Materials

Here's a collection of badly typed notes for the course, of which the naming conventions for the python file follows

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Current:
current, i, is measured as charge moved per time through an element
i = dq/dt  (dq=change in charge, dt=time interval)
current direction is direction of positve charge

1 A = 1 C/1 s


1 e = negative 1.602 x 10^-19 C of charge
1 C of charge = 6.24 x 10^18 e

-often see 'e' written as 'q'
-Charge is always a multiple of e charge
-Charge cannot be created or destroyed

Voltage:
amount of energy (J) to move 1 C of charge between points a and b
The voltage, or potential difference between points a and b, is the change in potential energy between 
these points divided by the charge that can be moved between these points.
potential difference Vab is potential energy difference divided by charge
Vab = Va - Vb = (Ua - Ub)/q
Here the subscript on U represents potential energy at that point, and q represents the charge moved between the points.


density of Si crystal is 5 x 10^22 Si atoms/cm^3
Doped Semiconductors:
ni - intrinsic concentration of e in conduction band
ni = 10^10 cm^-3 at room temp (T=300K) for Si
n - concentration of e in conduction band (measured in e/cm^3)
p - concentration of holes in valence band (measured in holes/cm^3)
np= ni^2
us this in doped SC but special cases

NsubA = acceptor impurtiy atom concentration (atoms/cm^3) p side
NsubD = donor impurtiy atom concentration (atoms/cm^3)    n side

"much greater than: >>" is at least 20 times greater 
case 1 : doped with only donor atoms to a concentration N sub D
    -if Nd >> ni, then n = Nd
    -if Nd is not >> ni, : n = ND/2 + [(ND)^2 /4 + ni^2]^(1/2)
once n is known, p is found from: p=ni^2 / n

case 2: doped with only acceptor atoms with concentration N sub A
    - if Na >> ni, then p = Na
    - if Na is not >> ni, then p = NA /2 + [(NA)^2 /4 + ni2]^(1/2)
once p is known, n is found from: n = ni^2 / p

DONOR ATOMS:
    Column/Group 5 - P (more common), As
    
ACCEPTOR ATOMS:
    Column/Group 3 - B (more common), Al

02.05a Dominant Carrier Type

n-type semiconductors: concentration of e in conduction band greater than holes in valence band. n>p (doped with donor atoms)
p-type semiconductors: concentration of holes in valence band greater than e in conduction band. p>n (dope with acceptor atoms)

SCs with both donor and acceptor are *compensated* SCs. can be either n or p type

-Majority/Minority Carriers-
    -In n-type, the majority carrier is an electron. The minority carrier is a hole.
    -In p-type, the majority carrier is a hole. The minority carrier is an electron.


Charge neutrality:  NsubD+ + p - NsubA- - n = 0

CASE 3: compensated SC, doped with donors to a level NsubD and acceptors to a level NsubA
a) if (NsubA - NsubD) >> nsubi the material is p-type. if NsubA >> NsubD, then p=NsubA and n=nsubi^2 / p;
    however if NsubA is NOT >> NsubD, p = (NsubA - NsubD)

b)if (NsubD - NsubA) >> nsubi the material is n-type. if NsubD >> NsubA, then n=NsubD and p=nsubi^2 / p;
    however if NsubD is NOT >> NsubA, p = (NsubD - NsubA)

c) if neither a nor b true, meaning |NsubD - NsubA| is NOT >> nsubi, then check which greater and apply equations...
    1) n-type, NsubD > NsubA and NsubD > nsubi
        n = ( (Nd - Na) + sqrt( (Nd-Na)^2 - 4ni^2) ) ) / 2, and p = ni^2 / n
    2) p-type, NsubA > NsubD and NsubA > nsubi
        p = ( (Na - Nd) + sqrt( (Na-Nd)^2 - 4ni^2) ) ) / 2, and n = ni^2 / n

d) if neither NsubA nor NsubD from c) are greater than nsubi, the material is dominated by the intrinsic carrier concentration and is not considered n or p type


02.05c Fermi Function and Fermi Level

-Fermi Function, a probability function, specifies under equilibrium conditions the probability that an available energy state at energy E will be occupied by an electron
    
    f(E) = 1 / (1 + e^( (E-EsubF) / kT) )

    where
    EsubF = Fermi Level
    E = energy level of interest
    k = Boltzman constant
    T = temp in K

Case 1. T = 0 K:


We refer to the energy level of instrinsic SCs as the intrinsic energy level, or Esubi


-E- Relationship between Esubi and EsubF

For n-type SC when NsubD >> NsubA and NsubD >> nsubi, the relationship between Esubi and EsubF is:
    
    EsubF - Esubi = kT * ln(NsubD/nsubi)

For p-type SC when NsubA >> NsubD and NsubA >> nsubi, the relationship between Esubi and EsubF is:
    
    EsubF - Esubi = kT * ln(NsubA/nsubi)

A pn-junction diode is designed to have dopant concentrations high enough to meet the conditions necessary for each of these equations to be valid.  
Therefore, we can apply each equation to the respective sides of the pn-junction to determine the placement of EF in each side.

Let's take an example of a pn-junction diode and calculate the positions of the Fermi energy levels. 
Suppose you are given a pn-junction diode with the p-side doped to a level of 10^16 /cm^3 and the n-side doped to a level of 10^17 /cm^3. 
 What is the position of the Fermi energy level relative to the intrinsic energy level on each side of the junction?  

On the p-side, we use the equation given above to calculate the position of EsubF and find:

    Esubi - EsubF = kT * ln(NsubA/nsubi) = 0.026 ln (10^16 / 10^10) = 0.359 eV

Notice that on the p-side, EsubF lies below Esubi, hence we take the difference Esubi - EsubF. and see that EsubF lies below Esubi by 0.359 eV

On the n-side, we find:
    
    EsubF - Esubi = kT * ln(NsubD/nsubi) = 0.026 ln (10^17 / 10^10) = 0.419 eV



02.08 Mobility, Resistivity, and Conductivity

Mobility - greek letter mu (|u), defined as a the slope of a line relating charge carrier (electron or hole) velocity in a material as a function of an 
    applied electric field (or voltage through a distance) across the material
UNITS: cm^2 / V-s
Impurity concentration N sub T. Usually just NsubA + NsubD, and zero if no impurities

Electron mobility, mu sub n : musubn = 52.2 + 1365/(1 + (NsubT/9.68*10^16)^0.68) cm^2 / V-s

Hole mobility, mu sub p: musubp = 44.9 + 426/(1 + (NsubT/2.23*10^16)^0.76)  cm^2 / V-s

Electrons can appear to move randomly with thermal motion in the lattice without an E-field. However, when an electric field is applied to the crystal,
the net electron motion is in a direction opposite to the electric field direction. This doesn't mean each electron moves in the same direction opposite 
to the electric field all the time; it simply means that the random motion of electrons in general results in a more directed motion. Electrons will still 
move in all directions, but the electric field causes more alignment to the random motion.


-Resistivity- 
A. Resistivity is a property of a material, represented by greek letter rho. SI unit Ohm-meter

    p (rho) = 1/ (  q(n*musubn + p*musubp))
    
B. Resistance. related to resistivity by geometry of material, cross section A, and length of material l

    R = p * l/A

Resistors - Types ( Fixed, variable(with pots or change under voltage or current conditions))
Materials - Carbon Film, Carbon Composition, Ceramic, Metal Film

-Conductivity-

Conductivity is the inverse of resistivity, greek symbol sigma                  sigma = 1 / p(rho)    unit: (Ohm-m)

conductivity can be calculated if we know e and hole mobilities, and e and hole carrier concentrations, using:
    sigma = q ( n*musubn + p*musubp)
    
Conductance (G) is the inverse of resistance, and has units of mhos (upside down omega), or Siemens (S). 
    G= 1/R

MODULE 3 - 

An electric field, E, is characterized by the energy required to move charge between locations in space, over a specific distance
E has units of energy / (charge-distance) = 1 J/(1 C-1m) = 1 V / 1 m. V/cm is often used
Charge carriers move, or drift, through a material under the influence of an applied electric field. The particles move with a velocity 
called the drift velocity, which depends on the strength of the electric field and the material the charge carrier is moving in.

The current resulting from the charge carrier motion in the presence of an electric field is called *drift current*. 
It’s the result of the movement of both holes and electrons.


++DRIFT CURRENT++

*Drift current density* is the current that moves through the entire volume of the material
represented by j

    j = Qv (C/cm^3)(cm/s) = A/cm^2
    
    A - Cross sectional area of material
    Q - charge density (Coulombs of charge in a unit volume). Can correspond to electron charge density(Qn), hole charge density(Qp), or total charge density
        Note that Qn = q * n    and Qp = q * p   and q = 1.602e-19
    v - velocity of charge in an electric field (carrier drift velocity)(either electrons or holes)

Remember that the symbol for current is i (or I).  When referring to current density, the symbol j (or J) is used.  Thus, i = j x A.

total *drift current density* is the sum of all hole and electron carrier drift current contributions since both charge carriers contribute to current,
and is given by:

    jsubtotal drift = jsubn drift + jsubp drift

drift velocity is v (or V sub d) = mu*E . So mobility is therefore mu = v/E

Mobility describes the ease of movement of carriers in an electric field.  The higher the mobility, the more easily a carrier can be transported through the material.
low electric field strength < 10^4 V/cm

Electron and hole drift currents are defined as:
    J sub n drift = Qsubn*vsubn   (electron drift current)
    J sub p drift = Qsubp*vsubp   (hole drift current)

We can rewrite these in terms of the relationship between drift velocity and electric field: 

    J sub n drift = Qnvn = (-qn )(-mun*E) = qn*mun*E   A/cm2
    
    J sub p drift = Qpvp = (+qp)(+mup*E) = qp*mup*E   A/cm2
    
Now the total drift current is the sum of the hole and electron drift currents:

    J sub Total drift = Jndrift + Jpdrift = q(n mun + p mup) E = sigma E   A/cm2

This relationship for drift current and applied electric field defines electrical conductivity:

    sigma = q(n mn + p mp)    (omega-cm)-1
    
And Resistivity p(rho) is the reciprical of conductivity:
    
    rho = 1/sigma = 1 / q(n mn + p mp) (omega-cm)

Calculation of conductivity (and therefore resistivity) can be simplified oftentimes when the sample is recognized to have a majority carrier much 
greater than the minority carrier concentration (and the intrinsic concentration)

Simplification 1: If the material is n-type and n>>p, then we can simplify conductivity and resistivity calculations to these equations:

    sigma = q(n*musubn) (omega-cm)-1
    
    rho = 1/sigma = 1/q(n*musubn) (omega-cm)
    
Simplification 2: If the material is p-type and p>>n, then we can simplify conductivity and resistivity calculations to these equations:

    sigma = q(p*musubp) (omega-cm)-1
    
    rho = 1/sigma = 1/q(p*musubp) (omega-cm)
    
And no matter which case applies, n-type, p-type, compensated, the full equation (from above) for resistivity can always be applied


--DIFFUSION--
In practical semiconductors, it is quite useful to create carrier concentration gradients, often by varying the dopant concentration (and/or the dopant type) across a region of semiconductor.

This gives rise to a diffusion current resulting from the natural tendency of carriers to move from high concentration regions to low concentration regions.  Diffusion current in a doped 
semiconductor is analogous to a gas moving across a room to evenly distribute itself across the volume. During diffusion, carriers move toward regions where they exist in a lower concentration, 
so diffusion current densities are proportional to the negative of the carrier gradient.

The diffusion current densities for holes, electrons, and total diffusion current density are given with the equations below:

    J sub n diffusion = (-q)*Dsubn*(-dn/dx) = qDsubn(dn/dx) A/cm^2
    
    J sub p diffusion = (q)*Dsubp*(-dp/dx) = -qDsubn(dn/dx) A/cm^2
    
    J sub Total diffusion = J sub n diffusion + J sub p diffusion = (-q)*Dsubn*(-dn/dx) + (q)*Dsubp*(-dp/dx) A/cm^2

where:
    -dn/dx and -dp/dx represent the concentration gradient of electrons and holes, respectively.
    
    Dsubn and Dsubp are the diffusivities, or diffusion constants, for electrons and holes in the semiconductor material, respectively
        The unit for diffusivity is cm^2 / s
    
--Einstein's Relationship: Mobility and Diffusivity --
Under conditions of equilibrium, when the net carrier currents must be zero (i.e. the sum of all currents present must be zero).  Note: this doesn't mean that each current contribution is zero, 
it just means that all current contributions sum to a total of zero current.  Assuming only the cases of drift and diffusion in a semiconductor, under equilibrium conditions, and focusing 
on the electrons, we can say that:

    J sub Total electron = J sub n drift + J sub n diffusion = q*n*musubn*E + q*Dsubn * (dn/dx) = 0
    
After some rearranging and applying some relationships outside the scope of this class, it can be found that:

    Dsubn/musubn = kT/q             Similarly for the case of holes:    Dsubp/musubp = kT/q

The quantity kT/q is used often in materials and electronic device fields and has a special name:  Thermal Voltage.  It is represented by VT, and at room temperature, it is estimated as VT = 26 mV


--GENERATION AND RECOMBINATION--

When a semiconductor is perturbed from the equilibrium state, an excess or deficit in the carrier concentration relative to the equilibrium values is created inside the semiconductor. 
Recombination-generation (R-G) is nature’s order-restoring mechanism, the means whereby the carrier excess or deficit inside the semiconductor is stabilized or eliminated.

Recombination - a process whereby electrons and holes are annihilated or destroyed.

Generation – a process whereby electrons and holes are created.

There are three main types of recombination and generation:

    Energy band-to-band
    R-G Center
    Auger recombination and it’s inverse, Impact Ionization generation.


+++++Module 4+++++
REMEMBER 
    NsubA = p side
    NsubD = n side


pn-junction at equilibrium
WsubD = depletion layer width, or space charge region width at equilibrium
Concentration of holes in SCR? 0
Concentration of electrons in SCR? 0
Concentration of negative charges in SCR? NsubA-
Concentration of positive charges in SCR? NsubD+


-Forward Bias-
The applied E (+V on the anode) is in the opposite direction of the built-in E, thus resulting in an E field across the SCR that is reduced to the difference between (E applied - E built-in)
When +V is high enough, the space charge region will be either significantly reduced or eliminated, allowing electrons to diffuse from n-side to p-side, and holes from p-side to n-side. 

-Reverse Bias-
The applied E (+V on cathode) adds to the build-in E, resulting in less diffusion of charge (insignificant current flow) and an increase of SCR size 



---04.05 pn-Junctions depletion regions, Energy Band diagrams---

-Voltage review-
Voltage is defined as the potential energy required to move charge through a distance, say from some point a to another point b

    1 V = 1 J / 1 C      That is, a potential difference of 1 V means it takes 1 J to move 1 C of charge through a distance from point a to b.

    electronvolt (eV) - defined as the energy needed to move 1 electron through a potential difference of 1 V

Example:
-If it takes 1 V potential difference to move 1 electron, how much energy in Joules does that correspond to?  We can now go back to our starting relationship that 1 V = 1 J/1C and calculate how many Joules corresponds to an energy of 1 eV.
-We start by noting that 1 electron has a charge of 1.602e-19 C.  Therefore, we can eliminate C from our equation by writing:
    
    1 V = (1 J/1 C) x (1.602x10-19 C/1 electron)  which simplifies to    1 V  = (1 J x 1.602e-19)/1 electron

-Now, bringing the 1 electron to the side with voltage, we rewrite this as:
    
    1 V x 1 electron = 1.602e-19 J         This defines the eV      1 eV = 1.602e-19 J

Calculating thermal energy, kT
   Thermal energy is defined as kT, where k is Boltzmann's constant and T is temperature in Kelvin.  
   At room temperature is kT is usually taken as 0.026 eV (because the definition of "room temperature" is subjective, we have standardized this to be the typical room temperature thermal energy).
   You'll want to memorize the thermal energy value at room temperature.  What is the room temperature thermal energy in Joules?
   
   1 eV = 1.602x10-19 J, so 0.026 eV x 1.602x10-19 J/eV = 4.16x10-21 J.

Calculating V sub T:
Example:
    If thermal energy, kT at room temperature is kT = 0.026 eV, what is the thermal voltage? We also know kT in J is 4.16e-21 J (from above ^^)
(stuff, do later)    
    
V sub T = 0.026 V


--Built-in Potential, V sub bi--

The built-in potential, V sub bi, is a measure of the shift in energy between the p-side and n-sides when equilibrium is established (usually measured in units of V)
    -It's directly related to the built-in electric field in the space charge region. This comes about due to a relationship between charge density and the electric field, given by Poisson's equation
    
Built-in Electric Field = -dV/dx (considering only one dimension)
    - dV/dx is simply the change in voltage with respect to position x

The built-in potential can be calculated if the doping concentrations of the p- and n-sides is known.  This is given by the relationship:

    V sub bi = V sub T ln((NsubA * NsubD)/nsubi^2)          (measured in V)
        where V sub T is the thermal energy

-Poisson's Equation-

d*xi/dx = rho(x)/(Ksubs * epsilon sub 0)            xi looks like big epsilon

rho(x) = q(p-n + (NsubD+) - (NsubA-))

xi = - dV/dx

dE/dx = - dV/dx         (E in eV)

    where xi is electric field (V/cm)
    Ksubs is dielectric contant of Si = 11.8
    epsilon sub 0 is permittivity of free space (8.85e-14 F/cm)
    rho(x) is charge density
    V is voltage
    E is energy (eV)


-Calculating depletion region widths-

We consider three widths associated with the depletion regions:  
    (1) the total depletion region width; 
    (2) the width of the depletion region on the p-side; and 
    (3) the width of the depletion region on the n-side
    
The depletion region widths on the p- and n-sides are determined by the doping concentration of their respective sides:  the more heavily doped a side is, the smaller the depletion region width on that side.      
*using x for these equations but it's really chi, I think*

p-side depletion width:
    x sub p = [(2*Ksubs * epsilon sub 0)/q  * NsubD/(NsubA(NsubA+NsubD))  *  Vsubbi]^1/2    

n-side depletion width:
    x sub n = [(2*Ksubs * epsilon sub 0)/q  * NsubA/(NsubD(NsubA+NsubD))  *  Vsubbi]^1/2
    
Total depletion width:
    W sub D = x sub p + x sub n =  [(2*Ksubs * epsilon sub 0)/q  * (NsubA+NsubD)/(NsubA*NsubD)  *  Vsubbi]^1/2

IF you're given an applied charge, subtract it from Vsubbi
    W sub D = x sub p + x sub n =  [(2*Ksubs * epsilon sub 0)/q  * (NsubA+NsubD)/(NsubA*NsubD)  *  (Vsubbi-VsubA)]^1/2

+++04.08++++

++++++++MODULE 5+++++++++++ PN JUNCTION DIODE CURRENT CALCULATIONS AND OPERATING REGIONS
5.03 Calculating current through pn junction diode

-Ideal diode current is calculated using:
    
    I = Isubs (e^(VsubA/VsubT) -1)
    
    Can be rewritten in solving for VsubA
    
    VsubA = VsubT * ln((I+Isubs)/Isubs)
    
    Or solving for Isubs
    
    Isubs = I / (e^(VsubA/VsubT) -1)
    
    where: 
    I - ideal diode current
    Isubs - reverse saturation current   (typically very small, less than 1e-12 A)
    VsubA - voltage applied
    VsubT - thermal voltage   (remember VsubT = kT/q)(assume 0.026 V for room temp) 
    
+++5.06+++
Energy band levels for 3 cases
(a) equilibrium - the Fermi energy level is flat (level) between the p- and n-sides;

(b) reverse bias - the n-side energy levels are moved down relative to the equilibrium position by an amount of energy corresponding to the energy of the applied voltage;

(c) forward bias - the n-side energy levels are moved up relative to the equilibrium position by an amount of energy corresponding to the energy of the applied voltage. 
    Under forward and reverse biasing, the Fermi energies no longer align between the p- and n-sides, but are shifted by the applied potential energy.
    

There is a relationship between the location of the Fermi energy level in a doped semiconductor and the intrinsic energy level.  

For an n-type semiconductor, when ND >> NA and ND >> ni, the relationship between EF and Ei is:
    
    EF - Ei = kT ln(ND/ni)          (units: eV)
    
For a p-type semiconductor, when NA >> ND and NA >> ni, the relationship between EF and Ei is:
    
    Ei - EF = kT ln(NA/ni)          (units: eV)
    
    where:
    EF - Fermi energy level
    Ei - intrinsic energy level
    
++++ MODULE 6+++++++

06.02 - Introduction to the MOSFET and MOS Capacitor

MOSFET - Metal Oxide Semiconductor Field Effect Transistor device 
In order to study the MOSFET, which is the main subject of Module 7, we need to first discuss one of the building block structures of the MOSFET, referred to as the MOS capacitor, or MOSCap for short

The MOSFET consists of a doped semiconductor (we'll used doped Si). 
This doped Si (referred to as the substrate, or body) will have a thin oxide layer on it, with a conductive layer (referred to as the gate) on top of the oxide layer.  
There are two regions on the substrate surface that are left uncovered by the oxide and gate.  
These two regions get doped with the opposite dopant type from that in the substrate. 
For example, if the substrate is doped p-type, these two regions are doped n-type.  These two regions create pn-junction diodes (see the MOSFET Structure figure below).  
During fabrication of the MOSFET, these two regions are placed on either side of the gate/oxide layers.  They are referred to as the Source and the Drain.  
Electrical contact is made to these pn-junction diodes on either side of the gate by placing an electrode directly on the top of those regions (referred to as Source and Drain electrodes).   
One final electrode is placed on the other side of the substrate (p-Si side in our example).  This is referred to as the Body electrode.  Overall, there are four electrodes:  Gate, Drain, Source, and Body. 



THE MOSCap

The MOS capacitor structure consists of a metal gate, an insulator (such as an oxide, SiO2), and a doped semiconductor.  

The function of the MOSCap layers are described in the figure below.  

    Gate conductor: easily moves carriers to the interface of the metal and the oxide layer causing charge buildup.

    Insulator:  keeps the charge at the gate separated from the semiconductor. An ideal insulator layer has no charge within it. 
        It is used to keep charge in the gate and the substrate separated; no current can flow through the insulator (since no charge can be within it).

    Doped semiconductor: allows charge carriers to move through it like we've seen when we studied drift and diffusion current densities in a doped semiconductor. 


+06.05+ Operating Regions of MOS Cap
When a voltage is applied across the MOSCap, it is applied across the gate to body.  
This is done by applying a voltage on vG with the body held at 0 V (or other fixed potential, as long as there is a potential difference between the gate and body).  
There are three distinct regions of operation of a MOSCap under applied potential, depending upon the range of the potential applied.  
These are:  accumulation, depletion, and inversion.  All of the regions of operation are defined based on how the applied gate voltage forces the substrate (body) charge carriers to respond.  

V sub g- applied voltage to gate

A) Accumulation
In the accumulation region, the gate voltage causes the majority carrier of the substrate to be attracted to the substrate/oxide interface 
(but remain on the side of the substrate since no carriers get into or move through the oxide).  

    V sub g < 0V (p-type Si)
    V sub g > 0V (n-type Si)
    

B) Depletion 

In the depletion region, a voltage is applied to the gate that will produce a space charge region, or depletion region, in the substrate near the substrate/oxide interface. 
The condition for operation in the depletion region is that the gate voltage is positive, but less than some 'threshold voltage', VT.
VT is defined as the potential at which the depletion region in the substrate has reached it's maximum width, and is no longer capable of forming more depletion region for charge balance.

VsubT > V sub g > 0V - p-Si substrate

C) Inversion
As we discussed for the case of Depletion operation, when vG is increased, the depletion region continues to grow until the p-Si can no longer grow a larger depletion width.  
It reaches a maximum, WDmax.  At this voltage (VT) where the depletion region is maximized, it is energetically more likely for the p-Si to break Si-Si bonds and generate electrons in the conduction band, than it is for more the depletion layer to grow.  
These free electrons generated are attracted towards the substrate/oxide interface to create charge balance with the gate.  
They create a thin layer of electrons at the interface.  This layer of electrons is referred to as an ‘inversion layer’.
Why is this line layer of electrons called an inversion layer? A: The layer is within a p-Si material. Since p-Si is dominated by a large hole concentration, a layer of electrons is considered an ‘inversion’ of the dominant carrier concentration.

V sub G > VsubT (p-type substrate)


+++++MODULE 7 ++++++++ Operating Regions of a MOSFET


The inversion channel in a p-type substrate MOSFET is comprised of electrons.  Therefore, we refer to p-type substrate MOSFETs as NMOS transistors
The 'N' helps us to remember the channel contains electrons; however it can be confusing because sometimes the 'N' is confused with the type of substrate.  Keep in mind, the 'N' refers to channel carriers!

Similarly, the n-type substrate MOSFETs are referred to as PMOS transistors.  The 'P' helps us remember the channel contains holes.

A. Overview of MOSFET Operation

The inversion layer in the MOSCap spans the length of the MOSCap structure and bridges the source and drain regions together, making contact with each.
 Voltages applied to the source and drain control movement of this inversion layer charge between the source and the drain.  
Traditionally, the drain is defined as the electrode to which the channel carriers move when a voltage is applied between the source and drain.  
For example, if the drain has a positive voltage on it (with respect to the source voltage) electrons would be attracted to the drain.  
Therefore, a positive drain voltage is used when the inversion layer consists of electrons (i.e. the substrate is p-type).  
If the drain has a negative voltage on it, holes would be attracted to the drain.  
Therefore a negative drain voltage is used when the inversion layer consists of hole carriers (i.e. the substrate is n-type).


The magnitude of the applied voltage difference between the source and drain, which is often just abbreviated by V sub DS, determines the current measured through the inversion channel.  In this way, the transistor is operated:  we control the current collected through the channel by applied voltage.  

     V sub DS - applied voltage difference between source and drain

There are two ways applied voltage can influence the current measured:

    1. Magnitude of the gate voltage which determines how much charge is in the inversion layer.
    2. Magnitude of the drain to source voltage difference which determines how much of the inversion layer charge is measured.

There are three operating regions for a MOSFET, and one other region that results in no operation (referred to as 'Cutoff').  The three operating regions are:

    - linear (also called triode)
    - pinch-off
    - saturation (aslo called active region)

++ LINEAR (or TRIODE) ++

The linear region corresponds to a small potential difference between the source and drain (i.e. V sub DS is small). 
VDS must stay within certain voltage boundaries, that differ depending upon the doping type of the substrate:
    1. p-type substrate (NMOS transistor)
        Triode conditions: 

            0 <= V sub DS <= (V sub GS - V sub TN)
            
            where:
                V sub GS - gate to source voltage
                V sub TN - threshold V for the MOSCap to go into inversion

            The condition that V sub DS must be less than or equal to VsubGS - VsubTN comes from the need to keep the MOSCap in inversion
            VsubGS is the gate to source voltage.  This is the same as the voltage across the MOSCap from gate to body when the source is at the same
            potential as the body of the MOSCap structure.

    2. n-type substrate ( PMOS transistor)
        Triode conditions:
            
            0 <= |VsubDS| <= |VsubGS - VsubTP|

            where:
                V sub TP is the threshold voltage for the MOSCap


++ Pinch-Off ++

In the pinch-off region, the MOSFET is transitioning from the linear to the saturation region. 
It is a point where the inversion channel begins to get constricted at the end near the drain and current measurement starts to get reduced from the channel.
It is not usually a region in which circuits are designed to operate.
In the pinch-off region, the voltage at the pinch-off point is defined as VsubDSAT, the voltage at which saturation begins, 
and is related to VsubGS and the threshold voltage by the following:

    NMOS: VsubGS - VsubDSAT = VsubTN
    
    PMOS: |VsubGS - VsubDSAT| = VsubTP


++ Saturation ++

In the saturation region, the voltage between the source and drain has very little effect on the measured channel current. 
In an ideal MOSFET, an increase in VsubDS would not result in an increase in channel current.  
However, in reality, an increase in VsubDS will be accompanied by very small increases in channel current.

For a MOSFET to operate in saturation, VsubDS must again stay within certain voltage boundaries, that differ depending upon the doping type of the substrate:

    1. p-type substrate ( NMOS transistor)
        Saturation conditions:              0 <= (VsubGS - VsubTN) <= VsubDS

The condition that VsubDS must be greater than or equal to VsubGS - VsubTN comes from the need to keep the MOSCap in inversion.  
VsubTN is the threshold voltage for the MOSCap to go into inversion. VsubGS is the gate to source voltage.  
This is the same as the voltage across the MOSCap from gate to body when the source is at the same potential as the body of the MOSCap structure.



    2. n-type substrate (PMOS transistor)
        Saturation conditions:              0 <= |VsubGS - VsubTN| <= |VsubDS|
        where: VsubTP is the threshold voltage


++ Cutoff ++

While not listed as an operational category, it should be noted that Cutoff is the region when the gate voltage has been reduced enough that the MOSCap structure goes out of inversion and into depletion.  
The MOSFET does not work when the MOSCap is in depletion; there are no free carriers in the channel between the source and drain in that case.



*IMAGE:REGIONS OF MOSFET OPERATION

The channel is represented by a solid pink color in the cartoon images above.  
The depletion region that exists in the MOSCap structure (remember: the inversion channel forms after the depletion width has been maximized in the MOSCap) is shown under the channel and the fixed acceptor ions are labeled as - with a circle around them. 
Also shown are the depletion regions of the reverse biased pn-junction diodes at the source/substrate and drain/substrate regions.

In pinch-off, the voltage applied to the drain creates such a large depletion region in the pn-junction diode between the drain and body that the depletion region impinges on the channel. 
The result is that the channel starts to get deformed and no longer consists of an even layer of charge under the oxide.  
It develops a 'pinched' region at the edge of the drain.  The current that can be collected is diminished due to the smaller contact between the channel and the drain electrode - now just allowing current to flow easily at the pinched point.  
Other carriers have to drift from the channel across the depletion region formed in order to reach the drain.  

A further increase in the drain voltage causes the depletion region of the pn-junction diode to grow even larger, thus impinging on the channel far into the channel length.  
In the figure above, the impingement is almost halfway through the channel length.  
Once past the pinch-off point, current no longer increases when the drain voltage is increased.  
This is because there is no direct contact between the channel carriers and the drain electrode like there is in the linear and pinch-off regions (although there is a only a very small contact point in pinch off).  
The current no longer increases with an increase in drain voltage. Instead all carriers collected at the drain are transported via drift across the depletion region by the electric field set up by the drain to source voltage.


__EXAMPLES__

Example 1.  Do we need to apply a +VD or –VD in order to move the channel carriers towards the drain in a p-type body MOSFET?

Solution.  Since the channel has electrons, we need to apply a +VsubD to the drain compared to the source, which is at VsubS, in order for the electrons to go towards the drain.  

Example 2. Given p-Si, Vs = 0V , VG > VTN, VB = 0V, VD = ?, what happens if VD < 0V?

Solution, If we look only at the diode formed between the drain and body, we can see that it is forward biased:

 0V -----|p-si   |   n+-Si|---- -VD
        <---Electron motion
This means we’d measure a large current at the drain due to the forward biased diode between the drain and the body, as well as the current that may be coming from the channel.  Therefore, we are not able to isolate and measure the current in the channel.


++++ 7.05 MOSFET Current-Voltage Curves +++

Family of curves - shows the different Drain-source current (mu A) vs Drain-source voltage (V) for different values of gate voltage. 
It shows 4 points on the curves:
    - Cutoff
    - Linear Region
    - Pinch off
    - Saturation

A. Cutoff Region
This is the region where the applied gate voltage is too low to place the MOSCap portion of the MOSFET into inversion.  Here the drain current, IsubD, is 0
This means the gate voltage, VsubGS (to keep MOSCap out of inversion) must be:  VsubGS ≤ VsubTN for NMOS or VsubGS ≥ VsubTP for PMOS.


B. Linear Region
In the low drain-to-source voltage region, the drain current displays a linear relationship to the drain-to-source voltage.  The equations and conditions depend on whether the MOSFET is an NMOS or PMOS.

NMOS, for (VsubGS - VsubTN) ≥ VsubDS ≥ 0:

    I sub D = K sub n (VsubGS - VsubTN - VsubDS/2) VsubDS

PMOS, for |VsubGS - VsubTP| ≥ |VsubDS| ≥ 0:

    I sub D = K sub p (VsubGS - VsubTN - VsubDS/2) VsubDS

In these equations, Kn and Kp are material dependent parameters, that are defined by the dimensions of the MOSFET channel (length and width), the mobility of the carrier in the channel, and the capacitance of the oxide layer.


C. Saturation Region
In the higher drain-to-source voltage region, the drain current displays a linear relationship to the drain-to-source voltage.  The equations and conditions depend on whether the MOSFET is an NMOS or PMOS.

NMOS, for VsubDS ≥ (VsubGS - VsubTN) ≥ 0:

    IsubD = Ksubn/2 * (VsubGS - VsubTN)^2 * (1 + lamba*VsubDS)


PMOS, for |VsubDS| ≥ |VsubGS - VsubTP| ≥ 0:

    IsubD = Ksubp/2 * (VsubGS - VsubTP)^2 * (1 + lamba*|VsubDS|)

In these equations, Kn and Kp are the same material dependent parameters as described above for the linear region equations.  
The parameter lambda is the channel-length modulation parameter that is empirically determined and depends on the channel length.  
Typical values for lambda are between 0 V^-1 ≤ l ≤ 0.2 V^-1.   


D. Threshold Voltage
The threshold voltage (VsubTN and VsubTP) for a given MOSFET depends on several factors, including any source to body potential difference (VsubSB, which is present when the source and body electrodes are not held at the same potential), the surface potential phi sub F (which represents the voltage across the substrate depletion layer at pinch-off), and a parameter, gamma, which is referred to as the body effect parameter, with units of V^1/2.  The equations for threshold voltage are given below.

NMOS:
    VsubTN = VsubTO + gamma*(sqrt(VsubSB + 2*phisubF)-sqrt(2*phisubF))

PMOS:
    VsubTP = VsubTO + gamma*(sqrt(VsubBS + 2*phisubF)-sqrt(2*phisubF))

where: 

    V sub TO - threshold voltage that is present when the source and body are tied to the same potential i.e. the potential difference between the body and source is 0 V)

In all of our discussions of the MOSFET and MOSCap threshold voltage, we have assumed that VTN = VTO (and VTP = VTO for the PMOS)
This has been the threshold voltage we've assumed so far. Notice that the source to body voltage in these threshold voltage equations is represented by VSB in the NMOS, but VBS in the PMOS devices.  
This notational difference is due to the sign convention associated with the electrode subscripts on a voltage and is something you will learn about in circuits and microelectronics.

The threshold voltage was also discussed in Module 07.02 C, and it is worth mentioning here that we can find the threshold voltage from the I-V curves by noting where pinch off occurs, VsubDSat, and defining the threshold voltages for NMOS and PMOS devices:

NMOS:  VsubGS - VsubDSat = VsubTN

PMOS:  |VsubGS - VsubDSat| = VsubTP

