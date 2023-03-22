#include "include/traj.h"
void generateFromFile(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    ifstream infile("input/pos.txt");
    if(!infile.is_open()){
        cerr << "Error reading file (sad)" << endl;
        exit(1);
    }
    string line;
    istringstream liness;
    float scale;
    getline(infile, line);
    liness = istringstream(line);
    liness >> scale;
    float x, y, z, v;
    float r = n_space * scale * a0; // occupy 80% of the space (change if needed)
    for(int p=0; p < 2; p++){
        int n = 0;
        getline(infile, line); // the header
        cout << "Reading " << line << endl;
        while(n < n_partd && getline(infile, line)){
            // input should be in this format: x, y, z, vx, vy, vz
            // bounds: -1 < x, y, z < 1; any number for vx, vy, vz, but be realistic.
            liness = istringstream(line);
            liness >> x;
            liness >> y;
            liness >> z;
            liness >> v;
            pos0x[p][n] = x * r;
            pos1x[p][n] = pos0x[p][n] + v * r;
            liness >> v;
            pos0y[p][n] = y * r;
            pos1y[p][n] = pos0y[p][n] + v * r;
            liness >> v;
            pos0z[p][n] = z * r;
            pos1z[p][n] = pos0z[p][n] + v * r;
            if(n == 0) cout << pos0x[p][n] << " " << pos1x[p][n] << " " << pos0y[p][n] << " " << pos1y[p][n] << " " 
            << pos0z[p][n] << " " << pos1z[p][n] << endl;
            n++;

            q[p][n] = qs[p];
            m[p][n] = mp[p];
            nt[p] += q[p][n];

            if(n == n_partd) cout << "SAMPLE: " << x << " " << y << " " << z << " " << v << endl;
        }
    }
}
void generateLineOfElectrons(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    int particles = 1024;
    float r = n_space * a0;
    float step = 2.0f / particles;
    for(int i=0; i < n_partd; i++){
        int part = i % particles; //FINISHME: hack method, should generate them sequentially instead of ABCABCABC
        pos0x[0][i] = (part * step - 1.0) * r * 0.8f; //-1.8e3
        pos1x[0][i] = pos0x[0][i] + 0.8f * r * step / 4000.0f;
        pos0y[0][i] = pos0z[0][i] = pos1y[0][i] = pos1z[0][i] = 0;
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;
        q[0][i] = qs[0]; q[1][i] = qs[1];
        m[0][i] = mp[0]; m[1][i] = mp[1];
        nt[0] += q[0][i]; nt[1] += q[1][i];
        //cout << pos0x[0][i] << " " << r * step / (float)100 << endl;
        //cout << pos0y[0][i] << " " << pos1y[0][i] << " " << pos0z[0][i] << " " << pos1z[0][i] << endl;
    }
}
void generate8ElectronCorners(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    float c = n_space * a0 * 0.8f;
    float xs[] = {c, c, c, c, -c, -c, -c, -c},
    ys[] = {c, c, -c, -c, c, c, -c, -c},
    zs[] = {c, -c, c, -c, c, -c, c, -c};
    for(int i=0; i < n_partd; ++i)
    {
        int id = i%8;
        pos0x[0][i] = pos1x[0][i] = xs[id];
        pos0y[0][i] = pos1y[0][i] = ys[id];
        pos0z[0][i] = pos1z[0][i] = zs[id];
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;
    }
}
void generateElectronAtLocation(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    int particles = 1024;
    float r = n_space * a0;
    float step = 2.0f / particles;
    for(int i=0; i < n_partd; i++){
        //pos0x[0][i] = pos0y[0][i] = pos0z[0][i] = pos1x[0][i] = pos1y[0][i] = pos1z[0][i] = -a0 / 2.0f;
        pos0x[0][i] = pos1x[0][i] = 0.f;
        pos1x[0][i] = -a0 * 0.005f;
        pos0y[0][i] = pos1y[0][i] = 0.f;
        pos0z[0][i] = pos1z[0][i] = a0 * -.25;//a0 * -30;
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;
        q[0][i] = qs[0]; //q[1][i] = qs[1];
        m[0][i] = mp[0]; //m[1][i] = mp[1];
        nt[0] += q[0][i]; //nt[1] += q[1][i];
    }
    /* backup code
        for(int i=0; i < n_partd; ++i){
        pos0x[0][i] = pos1x[0][i] = a0 * (n_space_divx - 1) / 2.0 - 1e-5;
        pos0y[0][i] = pos1y[0][i] = a0 * (n_space_divy - 1) / 2.0 - 1e-5;
        pos0z[0][i] = pos1z[0][i] = a0 * (n_space_divz - 1) / 2.0 - 1e-5;
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;

        q[0][i] = qs[0];
        m[0][i] = mp[0];
        nt[0] += q[0][i];
    }
    */
}
void generateSphere(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    //https://stackoverflow.com/a/16128461
    float half_cell = a0/2.f;
    int particles = n_partd;
    float dlong = pi * (3.0f - sqrt(5.0f)), dz = 2.0 / n_partd, lg = 0.0f, z = 1.0f - dz / 2.0f;
    for(int i=0; i < n_partd; ++i){
        float r = sqrt(1.0 - z * z) * r0;
        pos0x[0][i] = pos1x[0][i] = cos(lg) * r;// + half_cell;
        pos0y[0][i] = pos1y[0][i] = sin(lg) * r;// + half_cell;
        pos0z[0][i] = pos1z[0][i] = z * r0;// + half_cell;
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;
        z -= dz;
        lg += dlong;

        q[0][i] = qs[0];
        m[0][i] = mp[0];
        nt[0] += q[0][i];
    }
}
void generate2Electrons(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    float half_cell = a0/2;
    for(int i=0; i < n_partd; ++i){
        pos0x[0][i] = pos1x[0][i] = i < n_partd/2 ? half_cell : -half_cell;
        pos0y[0][i] = pos1y[0][i] = half_cell;
        pos0z[0][i] = pos1z[0][i] = half_cell;
        pos0x[1][i] = pos0y[1][i] = pos0z[1][i] = pos1x[1][i] = pos1y[1][i] = pos1z[1][i] = 0;

        q[0][i] = qs[0];
        m[0][i] = mp[0];
        nt[0] += qs[0];
    }
}
void generateParticles(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt){
    //generateLineOfElectrons(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
    //generate8ElectronCorners(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
    generateElectronAtLocation(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
    //generateSphere(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
    //generate2Electrons(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
}
void generateEmptyField(float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx]){
    //technically, we don't need to do this because 0 is just blank
    for (unsigned int i = 0; i < n_space_divx; i++)
    {
        //float x = ((float)i - (float)(n_space_divx) / 2.0) / ((float)(n_space_divx));
        for (unsigned int j = 0; j < n_space_divy; j++)
        {
            //float y = ((float)j - (float)(n_space_divy) / 2.0) / ((float)(n_space_divy));
            for (unsigned int k = 0; k < n_space_divz; k++)
            {
                //float z = ((float)k - (float)(n_space_divz) / 2.0) / ((float)(n_space_divz));
                //float r = sqrt(x * x + y * y + z * z); unused

                Ee[0][k][j][i] = 0; // 1000+i*100;
                Ee[1][k][j][i] = 0; // 2000+j*100;
                Ee[2][k][j][i] = 0; // 3000+k*100;

                Be[0][k][j][i] = 0; // z/sqrt(x*x+y*y);
                Be[1][k][j][i] = 0; // z/sqrt(x*x+y*y);
                Be[2][k][j][i] = 0; // 1*z*z+1;
            }
        }
    }
}
void generateStripedEField(float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx]){
    for (unsigned int k = 0; k < n_space_divz; k++) for (unsigned int j = 0; j < n_space_divy; j++) for (unsigned int i = 0; i < n_space_divx; i++){
        Ee[0][k][j][i] = i < 16 ? 1.f : 1e4;
        Ee[1][k][j][i] = 0.f;
        Ee[2][k][j][i] = 0.f;
        Be[0][k][j][i] = Be[1][k][j][i] = Be[2][k][j][i] = 0.f; 
    } 
}
void generateConstantBField(float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx]){
    for (unsigned int k = 0; k < n_space_divz; k++) for (unsigned int j = 0; j < n_space_divy; j++) for (unsigned int i = 0; i < n_space_divx; i++){
        Be[0][k][j][i] = 0.f;
        Be[1][k][j][i] = 0.0015f;
        Be[2][k][j][i] = 0.f;
    } 
}
void generateField(float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx]){
    generateEmptyField(Ee, Be);
    //generateStripedEField(Ee, Be);
    //generateConstantBField(Ee, Be);
}