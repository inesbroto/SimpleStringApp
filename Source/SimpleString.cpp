/*
  ==============================================================================

    SimpleString.cpp
    Created: 12 Feb 2021 1:10:03pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "SimpleString.h"

//==============================================================================
SimpleString::SimpleString (NamedValueSet& parameters, double k) : k (k)
{
//Cajón
    L_x = 0.75;
    L_y = 0.25;
    E_cj = 3.8e9;
    v = 0.3;
    thick = 0.005;
    D = E * pow(thick, 3) / (12 * (1 - pow(v, 2)));
    rho_cj = 1150;
    kappa = sqrt(D / (rho * thick));
    sigma0_cj = 1;
    sigma1_cj = 0.05;
    k_cj = 1.0 / sr;
    h_min = 2.1 * sqrt(k_cj * (sigma1_cj + sqrt(kappa * kappa + sigma1_cj * sigma1_cj)));


    N_x = ceil(L_x / h_min);
    N_y = ceil(L_y / h_min);
    h_cj = std::min(L_x / N_x, L_y / N_y);
    mu = kappa * k / (h_cj * h_cj);
    S = 2 * sigma1 * k / (h_cj * h_cj);



    for (int i = 0; i < 3; ++i)
        uStates_cj[i] = std::vector<std::vector<double>>(N_x, std::vector<double>(N_y, 0.0));

    uNext_cj = &uStates_cj[0];
    u_cj     = &uStates_cj[1];
    uPrev_cj = &uStates_cj[2];



    // An (N_x, N_y) x 3 'matrix' containing the state of the system at all time-step
    //uStates_cj = std::vector<std::vector<std::vector<float>>>(
    //    3, std::vector<std::vector<float>>(
    //        N_x + 1, std::vector<float>(N_y + 1, 0.0f)
    //    )
    //);
    //u_cj = std::vector<std::vector<float>> (N_x + 1, std::vector<float>(N_y + 1, 0.0f));
    //uPrev_cj = u_cj;
    //uNext_cj = u_cj;

    //excitation parameters
    //int exc_x, exc_y, exc_dev;

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
    A04 = sigma0 * k - 1 + 4 * S;
    A05 = -S;


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
    
    /*  Make u pointers point to the first index of the state vectors.
        To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
             - n = 0 is u^{n+1},
             - n = 1 is u^n, and
             - n = 2 is u^{n-1}.
        Also see calculateScheme()
     */
    
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
    g.strokePath(visualiseState (g, 100), PathStrokeType(2.0f));

}

Path SimpleString::visualiseState (Graphics& g, double visualScaling)
{
    // String-boundaries are in the vertical middle of the component
    double stringBoundaries = getHeight() / 2.0;
    
    // initialise path
    Path stringPath;
    
    // start path
    stringPath.startNewSubPath (0, -u[1][0] * visualScaling + stringBoundaries);
    
    double spacing = getWidth() / static_cast<double>(N);
    double x = spacing;
    
    for (int l = 1; l <= N; l++) // if you don't save the boundaries use l < N
    {
        // Needs to be -u, because a positive u would visually go down
        float newY = -u[1][l] * visualScaling + stringBoundaries;
        
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
    // Short references to reduce pointer dereferencing verbosity
    auto& u     = *u_cj;
    auto& uPrev = *uPrev_cj;
    auto& uNext = *uNext_cj;

        // Main update loop
    for (int l = 3; l < N_x - 2; ++l)
    {
        for (int m = 3; m < N_y - 2; ++m)
        {
            uNext[l][m] =
                A00 * u[l][m] +
                A01 * (u[l + 1][m] + u[l - 1][m] + u[l][m + 1] + u[l][m - 1]) +
                A02 * (u[l + 1][m + 1] + u[l - 1][m + 1] + u[l + 1][m - 1] + u[l - 1][m - 1]) +
                A03 * (u[l + 2][m] + u[l - 2][m] + u[l][m + 2] + u[l][m - 2]) +
                A04 * uPrev[l][m] +
                A05 * (uPrev[l + 1][m] + uPrev[l - 1][m] + uPrev[l][m + 1] + uPrev[l][m - 1]);
        }
    }

    // Mass-spring-collision terms
    for (int i = 0; i < numMasses; ++i)
    {
        int lc = massPos[i].first;
        int mc = massPos[i].second;

        double eta = x[i] - u[lc][mc];
        double phi = (eta > 0.0) ? (1.0 / (nu[i] + 1.0)) * K_col[i] * std::pow(eta, nu[i] + 1.0) : 0.0;
        double phiPrime = (eta > 0.0) ? K_col[i] * std::pow(eta, nu[i]) : 0.0;

        double xNew = 2.0 * x[i] - xPrev[i]
                    + (k * k / M) * (-K_mass[i] * (x[i] - x0[i]) - phiPrime - damping[i] * (x[i] - xPrev[i]) / k);

        uNext[lc][mc] += (k * k / (rho_cj * thick * h_cj * h_cj)) * phiPrime;

        xPrev[i] = x[i];
        x[i] = xNew;
    }
}




void SimpleString::updateStates_cajon()
{
    auto temp = uPrev_cj;
    uPrev_cj = u_cj;
    u_cj = uNext_cj;
    uNext_cj = temp;
    /*
    std::vector<std::vector<float>>* tmp = uPrev_cj;
    uPrev_cj = u_cj;
    u_cj = uNext_cj;
    uNext_cj = tmp;

    OR

    std::swap(uPrev_cj, u_cj);
    std::swap(u_cj, uNext_cj);
    */
}



/*
void SimpleString::calculateScheme_cajon()
{
    for (int l = 3; l < N_x - 2; ++l)
    {
        for (int m = 3; m < N_y - 2; ++m)
        {
            uNext_cj[l][m] = A00 * u_cj[l][m]
                           + A01 * (u_cj[l + 1][m] + u_cj[l - 1][m] + u_cj[l][m + 1] + u_cj[l][m - 1])
                           + A02 * (u_cj[l + 1][m + 1] + u_cj[l - 1][m + 1] + u_cj[l + 1][m - 1] + u_cj[l - 1][m - 1])
                           + A03 * (u_cj[l + 2][m] + u_cj[l - 2][m] + u_cj[l][m + 2] + u_cj[l][m - 2])
                           + A04 * uPrev_cj[l][m]
                           + A05 * (uPrev_cj[l + 1][m] + uPrev_cj[l - 1][m] + uPrev_cj[l][m + 1] + uPrev_cj[l][m - 1]);
        }
    }

    // Handle collisions
    for (int i = 0; i < numMasses; ++i)
    {
        int lc = massPos[i].first;
        int mc = massPos[i].second;

        double eta = x[i] - u_cj[lc][mc];
        double phi = (eta > 0) ? (1.0 / (nu[i] + 1)) * K_col[i] * pow(eta, nu[i] + 1) : 0.0;
        double phiPrime = (eta > 0) ? K_col[i] * pow(eta, nu[i]) : 0.0;

        double xNew = 2 * x[i] - xPrev[i] + (k * k / M) *
                      (-K_mass[i] * (x[i] - x0[i]) - phiPrime - damping[i] * (x[i] - xPrev[i]) / k);

        uNext_cj[lc][mc] += (k * k / (rho * thick * h_cj * h_cj)) * phiPrime;

        xPrev[i] = x[i];
        x[i] = xNew;
    }
}


void SimpleString::calculateScheme_cajon()
{
    for (int l = 3; l < N_x - 2; ++l)
    {
        for (int m = 3; m < N_y - 2; ++m)
        {
            uNext_cj[l][m] = (2 - 20 * mu * mu - 4 * S) * u_cj[l][m]
                + (8 * mu * mu + S) * (u_cj[l + 1][m] + u_cj[l - 1][m] + u_cj[l][m + 1] + u_cj[l][m - 1])
                - 2 * mu * mu * (u_cj[l + 1][m + 1] + u_cj[l - 1][m + 1] + u_cj[l + 1][m - 1] + u_cj[l - 1][m - 1])
                - mu * mu * (u_cj[l + 2][m] + u_cj[l - 2][m] + u_cj[l][m + 2] + u_cj[l][m - 2])
                + (sigma0 * k - 1 + 4 * S) * uPrev_cj[l][m]
                - S * (uPrev_cj[l + 1][m] + uPrev_cj[l - 1][m] + uPrev_cj[l][m + 1] + uPrev_cj[l][m - 1]);
        }
    }

    for (int i = 0; i < numMasses; ++i)
    {
        int lc = massPos[i].first;
        int mc = massPos[i].second;

        double eta = x[i] - u_cj[lc][mc];
        double phi = (eta > 0) ? (1.0 / (nu[i] + 1)) * K_col[i] * pow(eta, nu[i] + 1) : 0.0;
        double phiPrime = (eta > 0) ? K_col[i] * pow(eta, nu[i]) : 0.0;

        double xNew = 2 * x[i] - xPrev[i] + k * k / M *
                      (-K_mass[i] * (x[i] - x0[i]) - phiPrime - damping[i] * (x[i] - xPrev[i]) / k);

        uNext_cj[lc][mc] += k * k / (rho * thick * h_cj * h_cj) * phiPrime;

        xPrev[i] = x[i];
        x[i] = xNew;
    }
}
    

void SimpleString::updateStates_cajon()
{
    auto temp = uStates_cj[2];
    uStates_cj[2] = uStates_cj[1];
    uStates_cj[1] = uStates_cj[0];
    uStates_cj[0] = temp;

    u_cj = uStates_cj[1];
    uPrev_cj = uStates_cj[2];
    uNext_cj = uStates_cj[0];
}
*/







void SimpleString::calculateScheme()
{

    for (int l = 2; l < N-1; ++l) // clamped boundaries
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l + 1] + u[1][l - 1]) + B2 * (u[1][l + 2] + u[1][l - 2])
                + C0 * u[2][l] + C1 * (u[2][l + 1] + u[2][l - 1]);
    
    u[0][1] = Bss * u[1][1] + B1 * (u[1][2] + u[1][0]) + B2 * u[1][3]
            + C0 * u[2][1] + C1 * (u[2][2] + u[2][0]);
    u[0][N-1] = Bss * u[1][N-1] + B1 * (u[1][N] + u[1][N-2]) + B2 * (u[1][N-3])
            + C0 * u[2][N-1] + C1 * (u[2][N] + u[2][N-2]);


    
}

void SimpleString::updateStates()
{
    // Do a pointer-switch. MUCH quicker than copying two entire state vectors every time-step.
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
}

void SimpleString::excite()
{
    //// Arbitrary excitation function (raised cosine) ////
    
    // width (in grid points) of the excitation
    double width = 10;
    
    // make sure we're not going out of bounds at the left boundary
    int start = std::max (floor((N+1) * excitationLoc) - floor(width * 0.5), 1.0);

    for (int l = 0; l < width; ++l)
    {
        // make sure we're not going out of bounds at the right boundary (this does 'cut off' the raised cosine)
        if (l+start > (clamped ? N - 2 : N - 1))
            break;
        
        u[1][l+start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width-1.0)));
        u[2][l+start] += 0.5 * (1 - cos(2.0 * double_Pi * l / (width-1.0)));
    }
    // Disable the excitation flag to only excite once
    excitationFlag = false;
}

void SimpleString::mouseDown (const MouseEvent& e)
{
    // Get the excitation location as a ratio between the x-location of the mouse-click and the width of the app
    excitationLoc = e.x / static_cast<double> (getWidth());
    
    // Activate the excitation flag to be used by the MainComponent to excite the string
    excitationFlag = true;
}
