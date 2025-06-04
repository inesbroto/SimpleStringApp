#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>


#include <unistd.h>

using namespace std;

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h" // https://github.com/mackron/dr_libs/blob/master/dr_wav.h

void save_wav(const char* filename, const std::vector<float>& samples, unsigned int sample_rate) {
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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: ./synth norm_x norm_y\n";
        return 1;
    }



    float norm_x = atof(argv[1]);  // normalized x in [0, 1]
    float norm_y = atof(argv[2]);  // normalized y in [0, 1]
    //usleep(5000000);
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

    vector<vector<float>> u(N_x + 1, vector<float>(N_y + 1, 0.0));
    vector<vector<float>> uPrev = u;
    vector<vector<float>> uNext = u;

    int exc_x = floor(N_x * norm_x);
    int exc_y = floor(N_y * norm_y);
    int exc_dev = 2;
    auto hann = hann2D(exc_dev * 2 + 1);

    for (int i = -exc_dev; i <= exc_dev; ++i)
        for (int j = -exc_dev; j <= exc_dev; ++j)
            if (exc_x + i >= 0 && exc_x + i <= N_x && exc_y + j >= 0 && exc_y + j <= N_y)
                u[exc_x + i][exc_y + j] = hann[i + exc_dev][j + exc_dev];

    uPrev = u;

    const int duration = 1;
    const int numSamples = duration * sr;
    vector<float> out(numSamples, 0.0);


    int numMasses = 12;
    vector<float> x(numMasses, 0.0), xPrev(numMasses, 0.0);
    vector<float> x0(numMasses, 0.01);
    double M = 0.05;
    vector<float> K_mass(numMasses), K_col(numMasses, 1e8), nu(numMasses, 1.7), damping(numMasses, 0.01);

    for (int i = 0; i < numMasses; ++i)
        K_mass[i] = 800 + 100.0 * i / (numMasses - 1);

    vector<pair<int, int>> massPos = {
        {int(N_x * 0.4), int(N_y * 0.2)}, {int(N_x * 0.2), int(N_y * 0.2)},
        {int(N_x * 0.3), int(N_y * 0.2)}, {int(N_x * 0.6), int(N_y * 0.2)},
        {int(N_x * 0.8), int(N_y * 0.2)}, {int(N_x * 0.7), int(N_y * 0.2)},
        {int(N_x * 0.4), int(N_y * 0.4)}, {int(N_x * 0.2), int(N_y * 0.4)},
        {int(N_x * 0.3), int(N_y * 0.4)}, {int(N_x * 0.6), int(N_y * 0.4)},
        {int(N_x * 0.8), int(N_y * 0.4)}, {int(N_x * 0.7), int(N_y * 0.4)}
    };

    std::cout<<"Cajon grid size (x,y) = " <<N_x<<" "<<N_y<< '\n';
    std::cout<<"Normalized excitation coords (x,y) = " <<norm_x<<" "<<norm_y<< '\n';

    std::cout<<"Grid excitation coords (x,y) = " <<exc_x<<" "<<exc_y<< '\n';
    std::cout<<"Mass1 grid coords (x,y) = " <<massPos[0].first<<" "<<massPos[0].second<<'\n';


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

    std::string filename = "click_" + to_string(norm_x).substr(0,4) + "_" + to_string(norm_y).substr(0,4) + ".wav";
    save_wav(filename.c_str(), out, sr);
    system(("aplay " + filename).c_str());


    
    bool save_sound = false; // set to true if you want to keep the .wav generated sound
    // Delete the file
    if (not save_sound and std::remove(filename.c_str()) != 0) {
        perror("Error deleting file");
    } else {
        std::cout << "File deleted successfully\n";
    }

    return 0;
}
