#include <iostream>
#include <fstream>
#include <mpi.h>
#include <string>
#include <cstring>
#include <cmath>
using namespace std;

struct F_lims
{
    int bottom;
    int top;
};
/*
int calc_lims_beads(int n, int rank, int c_size) // n - chain length, c_size - comm_size
{
    int necessary_m_size = (n*n-n)/2;
    float size_per_rank = float(necessary_m_size) / float(c_size);
    int matr_num_el;
    int el_b, el_e;
    
    if (floor(c_size) == k) matr_num_el = int(size_per_rank);
    else {
        if (rank == 0) {
            el_b = 0; 
            el_e = necessary_m_size - int(size_per_rank)*(c_size-1);
        }
        else {
            el_b = necessary_m_size - int(size_per_rank)*(c_size-1) + int(size_per_rank)*rank;
            el_e = el_b + int(size_per_rank);
        }
    }
    return el_b, el_e;
}
*/
float calc_rc(float **x, float **y, float **z, int numberFrames, int numberBeads)
{
    float rc = 0;
    int f = numberFrames-1;
    float dx, dy, dz;
    for (unsigned int i = 1; i < numberBeads; i++)
    {
        dx = x[i][f] - x[i-1][f];
        dy = y[i][f] - y[i-1][f];
        dz = z[i][f] - z[i-1][f];
        rc += sqrt(dx*dx + dy*dy + dz*dz);
    }
    return rc/(numberBeads-1);
}

F_lims calc_lims_frames(int n, int rank, int k)
{
    F_lims mylims = { n / k * rank, n / k * (rank + 1)};
    return mylims;
}

int qsum (int m)
{
    int s=0;
    for (unsigned int i = 0; i < m; i++)
    {
        s += m;
    }
    return s;
}

void cm_aver (float **x, float **y, float **z, int **cm, int **cm_count, int frame_b, int frame_e, int numbeads, float rc, float **dm)
{
    float dx, dy, dz, r;
    for (unsigned int i = 0; i < numbeads; i++)
    {
        for (unsigned int j = 0; j < numbeads; j++)
        {
            for (unsigned int k = frame_b; k < frame_e; k++)
            {
                cm_count[i][j] += 1;
                dx = x[i][k] - x[j][k];
                dy = y[i][k] - y[j][k];
                dz = z[i][k] - z[j][k];
                r = sqrt(dx*dx + dy*dy + dz*dz);
                dm[i][j] += r;
                if (r < rc) cm[i][j] += 1;
            }
        }
    }
}

int main(int argc, char *argv[])
{    
    int rank, commsize, ierr;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&commsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int id;
    int numberBeads, numberFrames;
    string emp, line;
    fstream in;

    string mypath;
    if (rank==0) {
        cout << "Your input:\n";
        for (unsigned int i = 0; i < argc; i++) cout << argv[i] << ' ';
        cout << endl;
        if (argc < 2) {
            ierr = 1;
            cerr << "Number of arguments is zero, please check --help and/or your input." << endl;
            MPI_Finalize();
            return 0;  
        }
    }        

    for (unsigned int i = 0; i < argc; i++)
    {
        if (rank==0 && (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)) {
            cout << "Program to calculate contact and distance matrices using lammpstrj file. 
            Output .dat file will be in the same folder with suffix \'_cm.dat\' and \'_dm.dat\', respectively." << 
            "-p or --path defines a path to a file to be analyzed\n" << endl;
            ierr = 3;
            MPI_Finalize();
            return 0; 
        }
        if (strcmp(argv[i], "-p") == 0|| strcmp(argv[i], "--path") == 0)
        {
            mypath = argv[i+1];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    float **x, **y, **z;
    int bead_begin, bead_end;

    in.open(mypath, ios::in);    
    if (!in.is_open()) {
        cerr<<"Somthing went wrong with opening lmp file" << endl;
        MPI_Finalize();
        return 15;
    }
    numberFrames = 0;
    while (getline(in, line)) {
        if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
        ++numberFrames;
    }
    in.close();
    if (rank==0) cout << "1st reading is over by " << rank << endl;
    in.open(mypath, ios::in);
    getline(in, emp);
    while (!(emp.find("ITEM: NUMBER OF ATOMS") != std::string::npos))
    {
        getline(in, emp);
    }
    in >> numberBeads;
    if (rank == 0) cout << "number of frames " << numberFrames << '\n' << "number of beads " << numberBeads << endl;
    F_lims bot_top_frames;
    bot_top_frames = calc_lims_frames(numberFrames, rank, commsize);

    x = new float*[numberBeads];
    for (int i = 0; i < numberBeads; i++)
    {
        x[i] = new float[numberFrames];
    }
    y = new float*[numberBeads];
    for (int i = 0; i < numberBeads; i++)
    {
        y[i] = new float[numberFrames];
    }
    z = new float*[numberBeads];
    for (int i = 0; i < numberBeads; i++)
    {
        z[i] = new float[numberFrames];
    }
    
    int frame=-1;
    while (getline(in, emp)) {
        frame++;
        while (!(emp.find("ITEM: ATOMS id type xu yu zu") != std::string::npos)) {
            if (!getline(in, emp)) break;
        }
        for (int i = 0; i < numberBeads; i++)
        {
            in >> id;
            id--;
            in >> emp >> x[id][frame] >> y[id][frame] >> z[id][frame];
        }
    }
    in.close();
    if (rank==0) cout << "2nd reading is over by " << rank << endl;
    if (rank == 0) cout << numberBeads << ' ' << commsize << endl;

    float** dm = new float*[numberBeads];
    int** cm = new int*[numberBeads];
    int** cm_count = new int*[numberBeads];
    for (int i = 0; i < numberBeads; i++)
    {
        dm[i] = new float[numberBeads];
        cm[i] = new int[numberBeads];
        cm_count[i] = new int[numberBeads];
        for (unsigned int j = 0; j < numberBeads; j++)
        {
            dm[i][j] = 0;
            cm[i][j] = 0;
            cm_count[i][j] = 0;
        }
    }
    float rcontact = calc_rc(x, y, z, numberFrames, numberBeads);
    cm_aver(x, y, z, cm, cm_count, bot_top_frames.bottom, bot_top_frames.top, numberBeads, rcontact, dm); //x, y, z, cm, cm_count, numberFrames, bead_begin, bead_end);
    float** dm_total = new float*[numberBeads];
    int** cm_total = new int*[numberBeads];
    int** cm_total_count = new int*[numberBeads];
    for (int i = 0; i < numberBeads; i++)
    {
        dm_total[i] = new float[numberBeads];
        cm_total[i] = new int[numberBeads];
        cm_total_count[i] = new int[numberBeads];
        for (unsigned j = 0; j < numberBeads; j++)
        {
            dm_total[i][j] = 0;
            cm_total[i][j] = 0;
            cm_total_count[i][j] = 0;
        }        
    }
    for (unsigned int i = 0; i < numberBeads; i++)
    {
        MPI_Reduce(dm[i], dm_total[i], numberBeads, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(cm[i], cm_total[i], numberBeads, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(cm_count[i], cm_total_count[i], numberBeads, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    if (rank == 0) {
        // write to file distances matrix
        fstream fout;
        fout.open(mypath + "_dm.dat", ios::out);
        for (unsigned int i = 0; i < numberBeads; i++)
        {
            for (unsigned j = 0; j < numberBeads; j++)
            {
                if (j<i) {
                    fout << dm_total[j][i] / cm_total_count[j][i] << '\t';
                }
                else {
                    fout << dm_total[i][j] / cm_total_count[i][j] << '\t';
                }
            }
            fout << '\n';
        }
        fout.close();

        // write to file contact matrix
        fout.open(mypath + "_cm.dat", ios::out);
        for (unsigned int i = 0; i < numberBeads; i++)
        {
            for (unsigned j = 0; j < numberBeads; j++)
            {
                if (j<i) {
                    fout << cm_total[j][i] << '\t';
                }
                else {
                    fout << cm_total[i][j] << '\t';
                }
            }
            fout << '\n';
        }
        fout.close();
    }
    cout << rank << " is done!" << endl;
    MPI_Finalize();
    return 0;    
}
