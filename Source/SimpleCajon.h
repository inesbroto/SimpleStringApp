/*
  ==============================================================================

    SimpleCajon.h
    Created: 12 Feb 2021 1:10:03pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class SimpleCajon  : public juce::Component
{
public:
    SimpleCajon (NamedValueSet& parameters, double k_cj);
    ~SimpleCajon() override;

    void paint (juce::Graphics&) override;
    void resized() override;
    
    // function to draw the state of the string
    //Path visualiseState (Graphics& g, double visualScaling);
    Path visualiseState_cj (Graphics& g, double visualScaling);

    //void calculateScheme();
    //void updateStates();
    void calculateScheme_cajon();
    void updateStates_cajon();
    
    //return u at the current sample at a location given by the length ratio

    double getOutput_cj (double Lratio)
    {
         //return u[1][static_cast<int> (round(N * Lratio))];
        int idx = round(Lratio*N_y) * (N_x + 1) + round(Lratio*N_x);
        //DBG(idx);
        float val = u_pointer_cj[1][idx];
        //DBG("value");
        //DBG(val);
        return u_pointer_cj[1][static_cast<int> (idx)];

    }
    /*
    double getOutput (double Lratio)
    {
        int idx = round(N * Lratio);
        //DBG(idx);
        float val =u[1][idx];
        //DBG("value");
        //DBG(val);
        return u[1][static_cast<int> (round(N * Lratio))];
    }*/

    
    //void excite();
    void excite2D();

    void mouseDown (const MouseEvent& e) override;
    
    bool shouldExcite() { return excitationFlag; };
    
private:
    
/*
//Simple string app
    // Model parameters
    double L, rho, A, T, E, I, cSq, kappaSq, sigma0, sigma1, lambdaSq, muSq, h, k;
    // Number of intervals (N+1 is number of points including boundaries)
    int N;
    
    // An (N+1) x 3 'matrix' containing the state of the system at all time-steps
    std::vector<std::vector<double>> uStates;
    
    // vector of pointers that point to state vectors
    std::vector<double*> u;

    // Scheme variables
    //  - Adiv for u^{n+1} (that all terms get divided by)
    //  - B for u^n
    //  - C for u^{n-1}
    //  - S for precalculated sigma terms
    
   double Adiv, B0, Bss, B1, B2, C0, C1, S0, S1;
    
*/
//cajón parameters
    // Model parameters
    //const int sr = 44100;
    double L_x, L_y;
    float E_cj, v, thick, D, rho_cj, kappa, sigma0_cj, sigma1_cj, k_cj, h_min, h_cj, mu, S;
    int N_x, N_y;

    // An (N_x, N_y) x 3 'matrix' containing the state of the system at all time-step

    std::vector<std::vector<double>> uStates_cj;

    std::vector<double*> u_pointer_cj;

    //excitation parameters
    int exc_x, exc_y, exc_dev;

    // Mass-spring-collision setup
    int numMasses;
    std::vector<float> x, xPrev, x0, K_mass,K_col,nu,damping;
    double M;
    std::vector<std::pair<int, int>> massPos;

        // Plate update scheme coefficients
    double A00, A01, A02, A03, A04, A05;
    
    

    // flag to tell MainComponent whether to excite the scheme or not
    bool excitationFlag = false;
    
    // initialise location of excitation
    double excitationLoc = 0.5;
    
    bool clamped = true;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimpleCajon)
};
