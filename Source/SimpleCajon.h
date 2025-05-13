/*
  ==============================================================================

    SimpleString.h
    Created: 12 Feb 2021 1:10:03pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class SimpleString  : public juce::Component
{
public:
    SimpleString (NamedValueSet& parameters, double k);
    ~SimpleString() override;

    void paint (juce::Graphics&) override;
    void resized() override;
    
    // function to draw the state of the string
    Path visualiseState (Graphics& g, double visualScaling);

    void calculateScheme();
    void updateStates();
    
    //return u at the current sample at a location given by the length ratio

    double getOutput (double Lratio)
    {
        return u[1][static_cast<int> (round(N * Lratio))];
    }
    
    void excite();
    void mouseDown (const MouseEvent& e) override;
    
    bool shouldExcite() { return excitationFlag; };
    
private:
    
//Simple string app TO DO translation

    // vector of pointers that point to state vectors
    std::vector<double*> u;
    /* Scheme variables
        - Adiv for u^{n+1} (that all terms get divided by)
        - B for u^n
        - C for u^{n-1}
        - S for precalculated sigma terms
    */
   double Adiv, B0, Bss, B1, B2, C0, C1, S0, S1;
    

//cajón parameters
    // Model parameters
    const int sr=44100; // probably can be removed and get JUCE sr

    double L_x, L_y;
    double E, v, thick, D, rho, kappa, sigma0,sigma1, k, h_min, h, mu, S;
    double N_x, N_y;

    // An (N_x, N_y) x 3 'matrix' containing the state of the system at all time-step
    std::vector<std::vector<std::vector<float>>> uStates;
    //std::vector<double*> u;

    std::vector<std::vector<float>> uPrev;
    std::vector<std::vector<float>> u;
    std::vector<std::vector<float>> uNext;

    //excitation parameters
    int exc_x, exc_y, exc_dev;

    // Mass-spring-collision setup
    int numMasses;
    std::vector<float> x, xPrev, x0, K_mass,K_col,nu,damping;
    double M;
    std::vector<std::pair<int, int>> massPos;


    
    

    // flag to tell MainComponent whether to excite the scheme or not
    bool excitationFlag = false;
    
    // initialise location of excitation
    double excitationLoc = 0.5;
    
    bool clamped = true;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimpleString)
};
