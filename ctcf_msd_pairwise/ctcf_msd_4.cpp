#include <iostream>
#include <fstream>
#include <mpi.h>
#include <string>
#include <cstring>
#include <cmath>
using namespace std;

void msd_1d(float *arr, float *msd, int numberFrames)
{
    int *count = new int[numberFrames];
    for (int i = 0; i < numberFrames; i++)
    {
        for (int j = i + 1; j < numberFrames; j++)
        {
            msd[j - i] += (arr[i] - arr[j]) * (arr[i] - arr[j]);
            count[j - i] += 1;
        }
    }
    for (int i = 1; i < numberFrames; i++)
    {
        if (count[i] != 0)
            msd[i] /= count[i];
        else
            msd[i] = 0;
    }
    delete[] count;
    count = NULL;
}

int main(int argc, char *argv[])
{
    int rank, commsize, ierr;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int numberBeads, numberFrames, skip_rows;
    string emp, line;
    string mypath;
    fstream in;
    float *ctcf_dst, *ctcf_msd, *joint_msd;
    skip_rows = 0;

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
            cout << "Program to calculate pairwise MSD using dsts file.\n"
                 << "Output MSD file will be in the same folder with suffix \'_msd_ctcf.dat\'."
                 << "\nInputs:\n"
                 << "-p, --path    defines a path to a file to be analyzed\n"
                 << endl;
            ierr = 3;
            MPI_Finalize();
            return 3;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--path") == 0)
        {
            mypath = argv[i + 1];
        }
        if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--skip") == 0)
        {
            skip_rows = atoi(argv[i + 1]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        in.open(mypath, ios::in);
        if (!in.is_open())
        {
            cerr << "Somthing went wrong with opening dsts file" << endl;
            MPI_Finalize();
            return 15;
        }
        numberFrames = -1; // because of the header line
        while (getline(in, line))
        {
            ++numberFrames;
        }
        in.close();
        numberFrames -= skip_rows;
    }
    MPI_Bcast(&numberFrames, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ctcf_dst = new float[numberFrames];
    ctcf_msd = new float[numberFrames];
    for (int i = 0; i < numberFrames; i++)
    {
        ctcf_dst[i] = 0;
        ctcf_msd[i] = 0;
    }
    if (rank == 0)
        cout << "1st reading is over by " << rank << endl;
    in.open(mypath, ios::in);
    getline(in, emp);
    int line_number = 0;
    while (line_number < skip_rows)
    {
        getline(in, line);
        ++line_number;
    }
    line_number = 0;
    while (getline(in, emp))
    {
        line_number++;
        for (int i = 0; i < 8; i++)
        {
            if (rank == i)
                in >> ctcf_dst[line_number];
            else
                in >> emp;
        }
    }
    in.close();
    if (rank == 0)
        cout << "2nd reading is over by " << rank << '\t' << numberFrames << endl;

    msd_1d(ctcf_dst, ctcf_msd, numberFrames);
    delete[] ctcf_dst;
    ctcf_dst = NULL;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        joint_msd = new float[numberFrames];
        for (int i = 0; i < numberFrames; i++)
            joint_msd[i] = 0;
    }
    MPI_Reduce(ctcf_msd, joint_msd, numberFrames, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] ctcf_msd;
    ctcf_msd = NULL;
    if (rank == 0)
    {
        for (int i = 0; i < numberFrames; i++)
            joint_msd[i] /= commsize;
        ofstream out;
        out.open(mypath + "_msd_ctcf_rare.dat", ios::out);
        for (int i = 0; i < numberFrames; i++)
            out << joint_msd[i] << endl;
        out.close();
    }
    if (rank == 0)
    {
        delete[] joint_msd;
        joint_msd = NULL;
    }
    cout << rank << " is done!" << endl;
    MPI_Finalize();
    return 0;
}