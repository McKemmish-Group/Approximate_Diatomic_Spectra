# Function for calculating rotational, vibrational and electronic spectra
# Date: 31/7/2018

# Main Function Rovibronic
def E_evr(ID):
    # Variable parameters:
    # Number of each vibrational and rotational levels for the entered number of electronic levels
    max_v = 4
    max_J = 30
    max_E = countInput(ID)    
    # Boltzmann constant and room temperature
    kT = 206 #cm^-1    
    # Threshold/Tolerance for intensity
    Threshold = 10e-5
    
    # Call function to get variables
    Te, we, wexe, Be, De, alpha_e, beta_e, gamma_e, Total_Electronic_Spin, Symbol, Reflection_Symmetry = getVar(ID, max_E)
    print (Te, we, wexe, Be, De, alpha_e, beta_e, gamma_e, Total_Electronic_Spin, Symbol, Reflection_Symmetry)
    # Calculate Energy Levels
    EN = EnergyLevels(ID, max_E, max_v, max_J, Te, we, wexe, Be, De, alpha_e, beta_e, gamma_e)
    
    # Calculate Transition Energies & Intensities
    dE, IE = TransitionElec(ID, we, Be, max_E, max_v, max_J, EN, kT, Threshold, Total_Electronic_Spin, Symbol, Reflection_Symmetry)
    
    # Stem plot spectra of Transition Energy vs Intensity
    plotSpectra(ID, max_E, max_J, max_v, dE, IE)
    
    return
# End main function


# Function to count number of electronic states in input file
def countInput(ID):
    x = open('{0}input.txt'.format(ID), 'r')    # Open file for reading input
    f = x.readline()
    
    numElecStates = 0
    while f:
        numElecStates = numElecStates+1
        f = x.readline()
    
    x.close()
    return numElecStates

# Function to Grab Variables
def getVar(ID, max_E):
    
    x = open('{0}input.txt'.format(ID), 'r')    # Open file for reading input
    f = x.readline()
    
    Te = [0]*max_E
    we = [0]*max_E
    wexe = [0]*max_E
    Be = [0]*max_E
    De = [0]*max_E
    alpha_e = [0]*max_E
    beta_e = [0]*max_E
    gamma_e = [0]*max_E
    Total_Electronic_Spin = [0]*max_E
    Symbol = [""]*max_E
    Reflection_Symmetry = [""]*max_E
    
    # Read and Split input variables (Diatomic Constants!!!)
    i = 0
    while f:
        g = f.split()
        
        Te[i] = float(g[0])
        we[i] = float(g[1])
        wexe[i] = float(g[2])
        Be[i] = float(g[3])
        De[i] = float(g[4])
        alpha_e[i] = float(g[5])
        beta_e[i] = float(g[6])
        gamma_e[i] = float(g[7])
        Total_Electronic_Spin[i] = int(g[8])
        Symbol[i] = str(g[9])
        Reflection_Symmetry[i] = str(g[10])
        
        i = i+1
        f = x.readline()
    
    x.close()
    return Te, we, wexe, Be, De, alpha_e, beta_e, gamma_e, Total_Electronic_Spin, Symbol, Reflection_Symmetry


# Function to Calculate Energy Levels
def EnergyLevels(ID, max_E, max_v, max_J, Te, we, wexe, Be, De, alpha_e, beta_e, gamma_e):
    y = open('{0}outputLEVELS.txt'.format(ID), 'w')
    
    Total = max_J * max_v * max_E
    EN = [0] * Total
    
    # Calculate Energy Levels
    # Electronic Levels
    for E in range(0, max_E):
        # Vibrational Levels
        for v in range(0, max_v):
            # Rotational Levels
            for J in range(0, max_J):
                En = Te[E] + we[E]*(v+1/2) - wexe[E]*((v+1/2)**2) + Be[E]*J*(J+1) - De[E]*(J**2*(J+1)**2) - alpha_e[E]*J*(J+1)*(v+1/2) - beta_e[E]*J**2*(J+1)**2*(v+1/2) + gamma_e[E]*J*(J+1)*(v+1/2)**2
                y.write('{0}\t{1}\t{2}\t{3}\n'.format(E,v,J,En))
                
                # Store vector of En
                EN[E*max_v*max_J+v*max_J+J] = En

    y.close()
    return EN

# Function for Electronic Transitions
def TransitionElec(ID, we, Be, max_E, max_v, max_J, EN, kT, Threshold, Total_Electronic_Spin, Symbol, Reflection_Symmetry):
    import math # For use of exponential function in calculating intensities
    z = open('{0}outputTRANSITIONS.txt'.format(ID),'w')
    
    Tot_Jv = max_v*max_J
    Total = math.factorial(max_E)*math.factorial(max_v)*3*max_J # All combinations e.g E!*v!*3J (dJ = -1, 0, +1) NOTE: Typically excess size
    dE = [0]*Total
    IE = [0]*Total
    
    z.write('\nFrequency\t\tIntensity\t\tLower\t\t\tUpper\n\t\t\t\t\t\tE\tv\tJ\tE\tv\tJ\n')
    # Counter for vector storage with factorial nature of pure electronic transitions
    i = 0
    for E_Lower in range(0,max_E):
        for E_Upper in range(E_Lower,max_E):
            ElecFactor = SymmetryCondition(E_Upper, E_Lower, Total_Electronic_Spin, Symbol, Reflection_Symmetry)
            if (ElecFactor != 0):
                for v_Lower in range(0,max_v):
                    for v_Upper in range(v_Lower,max_v):
                        VTF = VertTransitionFactor(we, Be, E_Upper, E_Lower, v_Upper)
                        for J in range(0,max_J):
                            IEcheck = VTF*ElecFactor*(2*J+1)*math.exp(-EN[E_Lower*Tot_Jv+v_Lower*max_J+J]/kT)
                            if (IEcheck >= Threshold):
                                # dJ = -1 P
                                if (J > 1 and v_Upper != v_Lower):
                                    dE[i] = EN[E_Upper*Tot_Jv+v_Upper*max_J+(J-1)] - EN[E_Lower*Tot_Jv+v_Lower*max_J+J]
                                    IE[i] = IEcheck
                                    z.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(dE[i], IE[i], E_Lower, v_Lower, J, E_Upper, v_Upper, J-1))
                                    i = i + 1
                                # dJ = 0 Q
                                if (E_Upper != E_Lower or v_Upper != v_Lower):
                                    dE[i] = EN[E_Upper*Tot_Jv+v_Upper*max_J+J] - EN[E_Lower*Tot_Jv+v_Lower*max_J+J]
                                    IE[i] = IEcheck
                                    z.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(dE[i], IE[i], E_Lower, v_Lower, J, E_Upper, v_Upper, J))
                                    i = i + 1
                                # dJ = +1 R
                                if (J < max_J - 1):
                                    dE[i] = EN[E_Upper*Tot_Jv+v_Upper*max_J+(J+1)] - EN[E_Lower*Tot_Jv+v_Lower*max_J+J]
                                    IE[i] = IEcheck
                                    z.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(dE[i], IE[i], E_Lower, v_Lower, J, E_Upper, v_Upper, J+1))
                                    i = i + 1
                            else:
                                IE[i] = 0
    # print(i,Total) # NOTE: significantly less non-zero elements than total vector size
    
    # Create new vectors without zero elements for returning out of the function
    dEnew = [0]*i
    IEnew = [0]*i
    j = 0
    while (j < i):
        dEnew[j] = dE[j]
        IEnew[j] = IE[j]
        j = j + 1
    
    z.close()
    return dEnew, IEnew

# Function to Check the Symmetry Conditions for Vibronic Transitions
def SymmetryCondition(E_Upper, E_Lower, Total_Electronic_Spin, Symbol, Reflection_Symmetry):
    check = 0 # Default 0, change to 1 if passes symmetry conditions
    
    # Check dLambda = 0, ±1 and no change in Reflection Symmetry (Sigma Only)
    if (Symbol[E_Lower] == 'Sigma'):
        if (Symbol[E_Upper] == 'Sigma' and Reflection_Symmetry[E_Upper] == Reflection_Symmetry[E_Lower]):
            check = 1
        if (Symbol[E_Upper] == 'Pi'):
            check = 1
    if (Symbol[E_Lower] == 'Pi'):
        if (Symbol[E_Upper] == 'Sigma'):
            check = 1
        if (Symbol[E_Upper] == 'Pi'):
            check = 1
        if (Symbol[E_Upper] == 'Delta'):
            check = 1
    if (Symbol[E_Lower] == 'Delta'):
        if (Symbol[E_Upper] == 'Pi'):
            check = 1
        if (Symbol[E_Upper] == 'Delta'):
            check = 1
        if (Symbol[E_Upper] == 'Phi'):
            check = 1
    if (Symbol[E_Lower] == 'Phi'):
        if (Symbol[E_Upper] == 'Delta'):
            check = 1
        if (Symbol[E_Upper] == 'Phi'):
            check = 1
        if (Symbol[E_Upper] == 'Gamma'):
            check = 1
    if (Symbol[E_Lower] == 'Gamma'):
        if (Symbol[E_Upper] == 'Phi'):
            check = 1
        if (Symbol[E_Upper] == 'Gamma'):
            check = 1
    
    # Check dS = 0, if not change check back to 0
    if (Total_Electronic_Spin[E_Upper] != Total_Electronic_Spin[E_Lower]):
        check = 0
    
    return check

# Function to Estimate Vertical Transition
def VertTransitionFactor(we, Be, E_Upper, E_Lower, v_Upper):
    import math # For use of exponential function in calculating intensities
    
    v_vert = (we[E_Upper]/2)*(1/(math.sqrt(Be[E_Lower])) - 1/(math.sqrt(Be[E_Upper])))**2 - 1/2
    VTF = math.exp(-(v_Upper-v_vert)**2)
    
    return VTF

# Function to Plot Stem Spectra
def plotSpectra(ID, max_E, max_J, max_v, dE, IE):
    import matplotlib.pyplot as plt # For plotting spectra
    
    # Matplotlib - Graphing
    plt.title('{} Spectra'.format(ID))
    plt.xlabel('Transition Frequency (cm\u207B\u00B9)')
    plt.ylabel('Intensity')
    plt.stem(dE, IE, markerfmt=" ", basefmt=" ")
    plt.savefig('{0}outputSPECTRA_E{1}_v{2}_J{3}.png'.format(ID, max_E, max_v, max_J))
    plt.show()
    
    return

#Line to run function
#E_evr('C12O16')

#E_evr('Al27O16')

#E_evr('P31O16')