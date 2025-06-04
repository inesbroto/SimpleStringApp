/*
  ==============================================================================

    SimpleString.cpp
    Created: 12 Feb 2021 1:10:03pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "SimpleString.h"
#include <cmath>


//==============================================================================
SimpleString::SimpleString (NamedValueSet& parameters, double k_cj) : k_cj (k_cj)
{
//Cajón

    // Initialise member variables using the parameter set
    L_x = *parameters.getVarPointer ("L_x");
    L_y = *parameters.getVarPointer ("L_y");
    E_cj = *parameters.getVarPointer ("E_cj");
    v = *parameters.getVarPointer ("v");
    thick = *parameters.getVarPointer ("thick");
    rho_cj = *parameters.getVarPointer ("rho_cj");
    sigma0_cj = *parameters.getVarPointer ("sigma0_cj");
    sigma1_cj = *parameters.getVarPointer ("sigma1_cj");

    //L_x = 0.75;
    //L_y = 0.25;
    //E_cj = 3.8e9;
    //v = 0.3;
    //thick = 0.005;
    D = E_cj * pow(thick, 3) / (12 * (1 - pow(v, 2)));
    //rho_cj = 1150;
    kappa = std::sqrt(D / (rho_cj * thick));
    //sigma0_cj = 1;
    //sigma1_cj = 0.05;
    //k_cj = 1.0 / sr;

    h_min = 2.1 * sqrt(k_cj * (sigma1_cj + sqrt(kappa * kappa + sigma1_cj * sigma1_cj)));


    N_x = ceil(L_x / h_min);
    N_y = ceil(L_y / h_min);

    h_cj = std::min(L_x / N_x, L_y / N_y);
    mu = kappa * k_cj / (h_cj * h_cj);
    S = 2 * sigma1_cj * k_cj / (h_cj * h_cj);

    uStates_cj = std::vector<std::vector<double>> (3,
        std::vector<double>((N_x + 1)*(N_y + 1), 0));

    //std::vector<double*> u_pointer_cj; NOT DO iT. It makes the pointer null and unaccessable.

      //Make u pointers point to the first index of the state vectors.
        //To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
        //     - n = 0 is u^{n+1},
        //     - n = 1 is u^n, and
        //     - n = 2 is u^{n-1}.
        //Also see calculateScheme()
    
    
    // Initialise pointer vector
    u_pointer_cj.resize (3, nullptr);
    
    for (int i = 0; i < 3; ++i)
        u_pointer_cj[i] = &uStates_cj[i][0];


    // Mass-spring-collision setup
    numMasses=6;
    double M;

    x = std::vector<float> (numMasses, 0.0);
    xPrev = std::vector<float> (numMasses, 0.0);
    x0 = std::vector<float> (numMasses, 0.01);
    M = 0.05;
    K_mass = std::vector<float> (numMasses);
    K_col = std::vector<float> (numMasses, 1e8);
    nu = std::vector<float> (numMasses, 1.7);
    damping = std::vector<float> (numMasses, 0.01);


    for (int i = 0; i < numMasses; ++i)
        K_mass[i] = 800 + 100.0 * i / (numMasses - 1);

    massPos = {
        {int(N_x * 0.2), int(N_y * 0.4)}, {int(N_x * 0.2), int(N_y * 0.2)},
        {int(N_x * 0.2), int(N_y * 0.3)}, {int(N_x * 0.2), int(N_y * 0.6)},
        {int(N_x * 0.2), int(N_y * 0.8)}, {int(N_x * 0.2), int(N_y * 0.7)}
    };

        // Plate update scheme coefficients
    A00 = 2 - 20 * mu * mu - 4 * S;
    A01 = 8 * mu * mu + S;
    A02 = -2 * mu * mu;
    A03 = -mu * mu;
    A04 = sigma0_cj * k_cj - 1 + 4 * S;
    A05 = -S;
/*
// Simple String

    // Initialise member variables using the parameter set
    L = *parameters.getVarPointer ("L");
    rho = *parameters.getVarPointer ("rho");
    A = *parameters.getVarPointer ("A");
    T = *parameters.getVarPointer ("T");
    E = *parameters.getVarPointer ("E");
    I = *parameters.getVarPointer ("I");
    sigma0 = *parameters.getVarPointer ("sigma0");
    sigma1 = *parameters.getVarPointer ("sigma1");
    
    // Calculate wave speed (squared)
    cSq = T / (rho * A);
    
    // Calculate stiffness coefficient (squared)
    kappaSq = E * I / (rho * A);

    double stabilityTerm = cSq * k * k + 4.0 * sigma1 * k; // just easier to write down below
    
    h = sqrt (0.5 * (stabilityTerm + sqrt ((stabilityTerm * stabilityTerm) + 16.0 * kappaSq * k * k)));
    N = floor (L / h);
    h = L / N; // recalculate h
    
    lambdaSq = cSq * k * k / (h * h);
    muSq = kappaSq * k * k / (h * h * h * h);
    

    // Initialise vectors
    uStates = std::vector<std::vector<double>> (3,
                                        std::vector<double>(N+1, 0));
    
    //  Make u pointers point to the first index of the state vectors.
      //  To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
      //       - n = 0 is u^{n+1},
      //       - n = 1 is u^n, and
      //       - n = 2 is u^{n-1}.
      //  Also see calculateScheme()
    
    
    // Initialise pointer vector
    u.resize (3, nullptr);
    
    // Make set memory addresses to first index of the state vectors.
    for (int i = 0; i < 3; ++i)
        u[i] = &uStates[i][0];
    
    // Coefficients used for damping
    S0 = sigma0 * k;
    S1 = (2.0 * sigma1 * k) / (h * h);
    
    // Scheme coefficients
    B0 = 2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1; // u_l^n
    Bss = 2.0 - 2.0 * lambdaSq - 5.0 * muSq - 2.0 * S1;
    B1 = lambdaSq + 4.0 * muSq + S1;                   // u_{l+-1}^n
    B2 = -muSq;                                        // u_{l+-2}^n
    C0 = -1.0 + S0 + 2.0 * S1;                         // u_l^{n-1}
    C1 = -S1;                                          // u_{l+-1}^{n-1}
    
    Adiv = 1.0 / (1.0 + S0);                           // u_l^{n+1}
    
    // Divide by u_l^{n+1} term
    B0 *= Adiv;
    Bss *= Adiv;
    B1 *= Adiv;
    B2 *= Adiv;
    C0 *= Adiv;
    C1 *= Adiv;
*/
}

SimpleString::~SimpleString()
{
    
}

void SimpleString::paint (juce::Graphics& g)
{
    // clear the background
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    
    // choose your favourite colour
    g.setColour(Colours::cyan);
    
    // draw the state
    //g.strokePath(visualiseState (g, 100), PathStrokeType(2.0f));
    g.strokePath(visualiseState_cj (g, 100), PathStrokeType(2.0f));


}


Path SimpleString::visualiseState_cj (Graphics& g, double visualScaling)
{
    // String-boundaries are in the vertical middle of the component
    double stringBoundaries = getHeight() / 2.0;
    
    // initialise path
    Path stringPath;
    
    // start path
    int idx_y = floor(N_y*0.5);
    stringPath.startNewSubPath (0, -u_pointer_cj[1][idx_y*(N_x+1)] * visualScaling + stringBoundaries);
    
    double spacing = getWidth() / static_cast<double>(N_x);
    double x = spacing;
    
    for (int l = 1; l <= N_x; l++) // if you don't save the boundaries use l < N
    {
        // Needs to be -u, because a positive u would visually go down
        float newY = -u_pointer_cj[1][idx_y*(N_x+1)+l] * visualScaling + stringBoundaries;
        
        // if we get NAN values, make sure that we don't get an exception
        if (std::isnan(newY))
            newY = 0;
        
        stringPath.lineTo (x, newY);
        x += spacing;
    }
    // if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing

    return stringPath;
}

void SimpleString::resized()
{

}


void SimpleString::calculateScheme_cajon()

{
        // Main update loop
    for (int l = 3; l < N_x - 2; ++l)
    {
        for (int m = 3; m < N_y - 2; ++m)
        {
            int idx_x = m * (N_x + 1) + l;

            u_pointer_cj[0][idx_x] = 
                A00 * u_pointer_cj[1][idx_x] +
                A01 * (u_pointer_cj[1][idx_x + 1] + u_pointer_cj[1][idx_x - 1] + u_pointer_cj[1][(m+1) * (N_x + 1) + l] + u_pointer_cj[1][(m-1) * (N_x + 1) + l]) +
                A02 * (u_pointer_cj[1][(m+1) * (N_x + 1) + l+1] + u_pointer_cj[1][(m+1) * (N_x + 1) + l-1] + u_pointer_cj[1][(m-1) * (N_x + 1) + l+1] + u_pointer_cj[1][(m-1) * (N_x + 1) + l-1]) +
                A03 * (u_pointer_cj[1][idx_x + 2] + u_pointer_cj[1][idx_x - 2] + u_pointer_cj[1][(m+2) * (N_x + 1) + l] + u_pointer_cj[1][(m-2) * (N_x + 1) + l]) +
                A04 * u_pointer_cj[2][idx_x] +
                A05 * (u_pointer_cj[2][idx_x + 1] + u_pointer_cj[2][idx_x - 1] + u_pointer_cj[2][(m+1) * (N_x + 1) + l] + u_pointer_cj[2][(m-1) * (N_x + 1) + l]);
        }
    }

            


    // Mass-spring-collision terms
    //for (int i = 0; i < numMasses; ++i)
    //{
    //    int lc = massPos[i].first;
    //    int mc = massPos[i].second;
    //    double eta = x[i] - u[lc][mc];
    //    double phi = (eta > 0.0) ? (1.0 / (nu[i] + 1.0)) * K_col[i] * std::pow(eta, nu[i] + 1.0) : 0.0;
    //    double phiPrime = (eta > 0.0) ? K_col[i] * std::pow(eta, nu[i]) : 0.0;
    //    double xNew = 2.0 * x[i] - xPrev[i]
    //                + (k * k / M) * (-K_mass[i] * (x[i] - x0[i]) - phiPrime - damping[i] * (x[i] - xPrev[i]) / k);
    //    uStates_cj[0][lc][mc] += (k * k / (rho_cj * thick * h_cj * h_cj)) * phiPrime;
    //    xPrev[i] = x[i];
    //    x[i] = xNew;
    //}
}



void SimpleString::updateStates_cajon()
{
    double* uTmp = u_pointer_cj[2];
    u_pointer_cj[2] = u_pointer_cj[1];
    u_pointer_cj[1] = u_pointer_cj[0];
    u_pointer_cj[0] = uTmp;
}



void SimpleString::excite2D()
{
    //// 2D Raised Cosine Excitation ////
    DBG("Exciting!");
    double width = 10.0;  // width of excitation in grid points (assumed square for simplicity)
    double excitationX = 0.5;
    double excitationY = 0.5;

    int radius = static_cast<int>(width / 2.0);
    int centreX = static_cast<int>((N_x + 1) * excitationX); // excitationX: ∈ [0, 1]
    int centreY = static_cast<int>((N_y + 1) * excitationY); // excitationY: ∈ [0, 1]

    for (int dy = -radius; dy <= radius; ++dy)
    {
        int y = centreY + dy;
        if (y < 1 || y > N_y - 2) continue;

        for (int dx = -radius; dx <= radius; ++dx)
        {
            int x = centreX + dx;
            if (x < 1 || x > N_x - 2) continue;

            // radial distance
            double r = std::sqrt(dx * dx + dy * dy);
            if (r > radius) continue;


            const double  double_Pi  = MathConstants<double>::pi;
            // 2D raised cosine profile
            double raisedCos = 0.5 * (1.0 - std::cos(2.0 * double_Pi * (r / width)));

            int idx = y * (N_x + 1) + x;
            //DBG(idx);
            //DBG(raisedCos);
            u_pointer_cj[1][idx] += raisedCos;
            u_pointer_cj[2][idx] += raisedCos;
        }
    }

    excitationFlag = false;
}

void SimpleString::mouseDown (const MouseEvent& e)
{
    // Get the excitation location as a ratio between the x-location of the mouse-click and the width of the app
    excitationLoc = e.x / static_cast<double> (getWidth());
    
    // Activate the excitation flag to be used by the MainComponent to excite the string
    excitationFlag = true;
}
