#include <iostream>
#include <fstream>
#include <mpi.h>
#include <string>
#include <cstring>
#include <cmath>
using namespace std;

void nullify(float *arr, int size)
{
    for (unsigned int i = 0; i < size; i++)
    {
        arr[i] = 0;
    }
}

void calculate_pairwise_msd(float **x, float **y, float **z, float *pairwise_msd, float *msd_count, int numbeads, int numframes, int numbead_b, int numbead_e)
{
    float *pairwise_track;
    pairwise_track = new float[numframes];
    for (unsigned int i = numbead_b; i < numbead_e; i++)
    {
        for (unsigned int i2 = i + 1; i2 < numbeads; i2++)
        {
            nullify(pairwise_track, numframes);
            for (unsigned int j = 0; j < numframes; j++)
            {
                pairwise_track[j] = sqrt(
                    (x[i][j] - x[i2][j]) * (x[i][j] - x[i2][j]) +
                    (y[i][j] - y[i2][j]) * (y[i][j] - y[i2][j]) +
                    (z[i][j] - z[i2][j]) * (z[i][j] - z[i2][j]));
            }
            for (unsigned int j1 = 0; j1 < numframes; j1++)
            {
                for (unsigned int j2 = j1; j2 < numframes; j2++)
                {
                    pairwise_msd[j2 - j1] += (pairwise_track[j2] - pairwise_track[j1]) * (pairwise_track[j2] - pairwise_track[j1]);
                    msd_count[j2 - j1] += 1;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, commsize, ierr;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int id, reg_size, unreg_size;
    int numberBeads, numberFrames;
    string emp, line;
    fstream in;

    bool lmptrj = false;
    bool xyztrj = false;
    bool pairwise_flag = true;
    string mypath;
    if (rank == 0)
    {
        cout << "Your input:\n";
        for (unsigned int i = 0; i < argc; i++)
            cout << argv[i] << ' ';
        cout << endl;
        if (argc < 2)
        {
            ierr = 1;
            cerr << "Number of arguments is zero, please check help and/or your input." << endl;
            MPI_Finalize();
            return 1;
        }
    }

    for (unsigned int i = 0; i < argc; i++)
    {
        if (rank == 0 && (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0))
        {
            cout << "Program to calculate MSD using lammpstrj file.\n"
                 << "Output MSD file will be in the same folder with suffix \'_msd_a.dat\'."
                 << "\nIt has two modes defining by -m or --mode key:\n"
                 << "-p, --path    defines a path to a file to be analyzed\n"
                 << "-lmp            takes lammpstrj format for input data " << endl;
            ierr = 3;
            MPI_Finalize();
            return 3;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--path") == 0)
        {
            mypath = argv[i + 1];
        }
        if (strcmp(argv[i], "-lmp") == 0)
            lmptrj = true;
        if (strcmp(argv[i], "-xyz") == 0)
            xyztrj = true;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    float **x, **y, **z;
    int bead_begin, bead_end;

    if (lmptrj)
    {
        in.open(mypath, ios::in);
        if (!in.is_open())
        {
            cerr << "Somthing went wrong with opening lmp file" << endl;
            MPI_Finalize();
            return 15;
        }
        numberFrames = 0;
        while (getline(in, line))
        {
            if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
                ++numberFrames;
        }
        in.close();
        if (rank == 0)
            cout << "1st reading is over by " << rank << endl;
        in.open(mypath, ios::in);
        getline(in, emp);
        while (!(emp.find("ITEM: NUMBER OF ATOMS") != std::string::npos))
        {
            getline(in, emp);
        }
        in >> numberBeads;
        if (rank == 0)
            cout << "number of frames " << numberFrames << '\n'
                 << "number of beads " << numberBeads << endl;
        if (numberBeads % commsize == 0)
        {
            reg_size = int(numberBeads / commsize);
            unreg_size = reg_size;
        }
        else
        {
            reg_size = int(numberBeads / commsize);
            unreg_size = numberBeads - (commsize - 1) * reg_size;
        }
        if (rank == 0)
        {
            bead_begin = 0;
            bead_end = unreg_size;
        }
        else
        {
            bead_begin = unreg_size + (rank - 1) * reg_size;
            bead_end = unreg_size + rank * reg_size;
        }

        x = new float *[numberBeads];
        for (int i = 0; i < numberBeads; i++)
        {
            x[i] = new float[numberFrames];
        }
        y = new float *[numberBeads];
        for (int i = 0; i < numberBeads; i++)
        {
            y[i] = new float[numberFrames];
        }
        z = new float *[numberBeads];
        for (int i = 0; i < numberBeads; i++)
        {
            z[i] = new float[numberFrames];
        }

        int frame = -1;
        while (getline(in, emp))
        {
            frame++;
            while (!(emp.find("ITEM: ATOMS id type xu yu zu") != std::string::npos))
            {
                if (!getline(in, emp))
                    break;
            }
            for (int i = 0; i < numberBeads; i++)
            {
                in >> id;
                id--;
                if (id >= bead_begin && id < bead_end)
                {
                    in >> emp >> x[id][frame] >> y[id][frame] >> z[id][frame];
                }
                else
                {
                    in >> emp >> emp >> emp >> emp;
                }
            }
        }
        in.close();
        if (rank == 0)
            cout << "2nd reading is over by " << rank << endl;
    }

    if (rank == 0)
        cout << numberBeads << ' ' << reg_size << ' ' << unreg_size << ' ' << commsize << endl;
    if (pairwise_flag)
    {
        float *msd_pairwise = new float[numberFrames];
        float *msd_pairwise_count = new float[numberFrames];
        nullify(msd_pairwise, numberFrames);
        nullify(msd_pairwise_count, numberFrames);
        calculate_pairwise_msd(x, y, z, msd_pairwise, msd_pairwise_count, numberBeads, numberFrames, bead_begin, bead_end);

        float *msd_total = new float[numberFrames];
        float *msd_total_count = new float[numberFrames];
        nullify(msd_total, numberFrames);
        nullify(msd_total_count, numberFrames);
        MPI_Reduce(msd_pairwise, msd_total, numberFrames, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(msd_pairwise_count, msd_total_count, numberFrames, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            fstream fout;
            fout.open(mypath + "_msd_p.dat", ios::out);
            for (int i = 1; i < numberFrames; i++)
            {
                if (msd_total_count[i] != 0)
                    fout << i << '\t' << msd_total[i] / msd_total_count[i] << endl;
                else
                    fout << i << '\t' << '0' << endl;
            }
            fout.close();
        }
    }

    cout << rank << " is done!" << endl;
    MPI_Finalize();
    return 0;
}
