#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;





#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h" // Download from https://github.com/mackron/dr_libs/blob/master/dr_wav.h

void save_wav(const char* filename, const std::vector<float>& samples, unsigned int sample_rate)
{
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = 1;
    format.sampleRate = sample_rate;
    format.bitsPerSample = 32;

    drwav wav;
    if (drwav_init_file_write(&wav, filename, &format, NULL)) {
        drwav_write_pcm_frames(&wav, samples.size(), samples.data());
        drwav_uninit(&wav);
    }
}




// Hanning window function
vector<vector<double>> hann2D(int size) {
    vector<double> hann1D(size);
    for (int i = 0; i < size; ++i)
        hann1D[i] = 0.5 * (1 - cos(2 * M_PI * i / (size - 1)));

    vector<vector<double>> hann(size, vector<double>(size));
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            hann[i][j] = hann1D[i] * hann1D[j];
    return hann;
}






int main() {
    const int sr = 44100;
    double L_x = 0.75, L_y = 0.25;
    double E = 3.8e9, v = 0.3, thick = 0.005;
    double D = E * pow(thick, 3) / (12 * (1 - pow(v, 2)));
    double rho = 1150;
    double kappa = sqrt(D / (rho * thick));
    double sigma0 = 1, sigma1 = 0.05;
    double k = 1.0 / sr;
    double h_min = 2.1 * sqrt(k * (sigma1 + sqrt(kappa * kappa + sigma1 * sigma1)));

    int N_x = ceil(L_x / h_min);
    int N_y = ceil(L_y / h_min);
    double h = min(L_x / N_x, L_y / N_y);
    double mu = kappa * k / (h * h);
    double S = 2 * sigma1 * k / (h * h);


    std::cout<<"L_x "<<L_x<<"L_y "<<L_y<<'\n';;//0.75, L_y;//0.25;
    std::cout<<E<<" v: "<<v<<" thick :"<<thick<<'\n';//3.8e9, v<<'\n';//0.3, thick<<'\n';//0.005<<'\n';
    std::cout<<D<<'\n';//E * pow(thick, 3) / (12 * (1 - pow(v, 2)))<<'\n';
    std::cout<<rho<<'\n';//1150<<'\n';
    std::cout<<kappa<<'\n';//sqrt(D / (rho * thick))<<'\n';
    std::cout<<" sigma0 "<<sigma0<<" sigma1 "<<sigma1<<'\n';//1, sigma1<<'\n';//0.05<<'\n';
    std::cout<<k<<'\n';//1.0 / sr<<'\n';
    std::cout<<"hmin "<<h_min<<'\n';//2.1 * sqrt(k * (sigma1 + sqrt(kappa * kappa + sigma1 * sigma1)))<<'\n';
    std::cout<<N_x<<'\n';//ceil(L_x / h_min)<<'\n';
    std::cout<<N_y<<'\n';//ceil(L_y / h_min)<<'\n';
    std::cout<<h<<'\n';//min(L_x / N_x, L_y / N_y)<<'\n';
    std::cout<<"mu "<<mu<<'\n';//kappa * k / (h * h)<<'\n';
    std::cout<<"S "<<S<<'\n';//2 * sigma1 * k / (h * h);


    vector<vector<float>> u(N_x + 1, vector<float>(N_y + 1, 0.0));
    vector<vector<float>> uPrev = u;
    vector<vector<float>> uNext = u;

    int exc_x = floor(N_x * 0.6);
    int exc_y = floor(N_y * 0.6);
    int exc_dev = 2;
    auto hann = hann2D(exc_dev * 2 + 1);

    for (int i = -exc_dev; i <= exc_dev; ++i)
        for (int j = -exc_dev; j <= exc_dev; ++j)
            u[exc_x + i][exc_y + j] = hann[i + exc_dev][j + exc_dev];

    uPrev = u;

    const int duration = 1;
    const int numSamples = duration * sr;
    vector<float> out(numSamples, 0.0);

    // Mass-spring-collision setup
    int numMasses = 6;
    vector<float> x(numMasses, 0.0), xPrev(numMasses, 0.0);
    vector<float> x0(numMasses, 0.01);
    double M = 0.05;
    vector<float> K_mass(numMasses), K_col(numMasses, 1e8), nu(numMasses, 1.7), damping(numMasses, 0.01);

    for (int i = 0; i < numMasses; ++i)
        K_mass[i] = 800 + 100.0 * i / (numMasses - 1);

    vector<pair<int, int>> massPos = {
        {int(N_x * 0.2), int(N_y * 0.4)}, {int(N_x * 0.2), int(N_y * 0.2)},
        {int(N_x * 0.2), int(N_y * 0.3)}, {int(N_x * 0.2), int(N_y * 0.6)},
        {int(N_x * 0.2), int(N_y * 0.8)}, {int(N_x * 0.2), int(N_y * 0.7)}
    };

    for (int n = 0; n < numSamples; ++n) {
        for (int l = 3; l < N_x - 2; ++l) {
            for (int m = 3; m < N_y - 2; ++m) {
                uNext[l][m] = (2 - 20 * mu * mu - 4 * S) * u[l][m]
                    + (8 * mu * mu + S) * (u[l + 1][m] + u[l - 1][m] + u[l][m + 1] + u[l][m - 1])
                    - 2 * mu * mu * (u[l + 1][m + 1] + u[l - 1][m + 1] + u[l + 1][m - 1] + u[l - 1][m - 1])
                    - mu * mu * (u[l + 2][m] + u[l - 2][m] + u[l][m + 2] + u[l][m - 2])
                    + (sigma0 * k - 1 + 4 * S) * uPrev[l][m]
                    - S * (uPrev[l + 1][m] + uPrev[l - 1][m] + uPrev[l][m + 1] + uPrev[l][m - 1]);
            }
        }

        for (int i = 0; i < numMasses; ++i) {
            int lc = massPos[i].first;
            int mc = massPos[i].second;

            double eta = x[i] - u[lc][mc];
            double phi = (eta > 0) ? (1.0 / (nu[i] + 1)) * K_col[i] * pow(eta, nu[i] + 1) : 0.0;
            double phiPrime = (eta > 0) ? K_col[i] * pow(eta, nu[i]) : 0.0;

            double xNew = 2 * x[i] - xPrev[i] + k * k / M * (-K_mass[i] * (x[i] - x0[i]) - phiPrime - damping[i] * (x[i] - xPrev[i]) / k);
            uNext[lc][mc] += k * k / (rho * thick * h * h) * phiPrime;
            xPrev[i] = x[i];
            x[i] = xNew;
        }

        uPrev = u;
        u = uNext;

        out[n] = u[int(N_x * 0.7)][int(N_y * 0.6)];
    }
    //system("aplay tom.wav");

    /*
    // Save output to file
    ofstream fout("plate_output.txt");
    for (double sample : out)
        fout << sample << "\n";
    fout.close();

    save_wav("tom_center.wav", out, sr);
    system("aplay tom.wav");
    cout << "Simulation complete. Output saved to 'plate_output.txt'.\n";
    return 0;
    */
}

