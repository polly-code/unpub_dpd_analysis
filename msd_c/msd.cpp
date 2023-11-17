#include <iostream>
#include <fstream>
#include <mpi.h>
#include <string>
#include <cstring>
using namespace std;

void msd_aver(float **x, float **y, float **z, float *msd, float *msd_count, int numframes, int numbead_b, int numbead_e)
{
    float dx, dy, dz;
    for (int i = numbead_b; i < numbead_e; i++)
    {
        for (int j = 0; j < numframes; j++)
        {
            for (int k = j + 1; k < numframes; k++)
            {
                dx = x[i][k] - x[i][j];
                dy = y[i][k] - y[i][j];
                dz = z[i][k] - z[i][j];
                msd[k - j] += dx * dx + dy * dy + dz * dz;
                msd_count[k - j] += 1;
            }
        }
    }
}

void msd_single(float **x, float **y, float **z, float **msd, float **msd_count, int numframes, int numbead_b, int numbead_e)
{
    float dx, dy, dz;
    for (int i = numbead_b; i < numbead_e; i++)
    {
        for (int j = 0; j < numframes; j++)
        {
            for (int k = j + 1; k < numframes; k++)
            {
                dx = x[i][k] - x[i][j];
                dy = y[i][k] - y[i][j];
                dz = z[i][k] - z[i][j];
                msd[i][k - j] += dx * dx + dy * dy + dz * dz;
                msd_count[i][k - j] += 1;
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
    int start_frame = 0;
    string emp, line;
    fstream input_file;
    bool aver = false;
    bool single = false;
    bool lmptrj = false;
    bool xyztrj = false;
    bool more_cols_flag = false;
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
                 << "Output MSD file will be in the same folder with suffix \'_msd_[a/s].dat\'."
                 << "\nIt has two modes defining by -m or --mode key:\n"
                 << "\taver - which is average MSD over all beads,\n"
                 << "\tsingle - MSD per every single bead\n"
                 << "-p, --path     defines a path to a file to be analyzed\n"
                 << "-s, --start    start from this frame\n"
                 << "-mc,           enables parser of additional 3 columns"
                 << "-xyz           takes xyz format for input data\n"
                 << "-lmp           takes lammpstrj format for input data " << endl;
            ierr = 3;
            MPI_Finalize();
            return 3;
        }
        if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mode") == 0)
        {
            if (strcmp(argv[i + 1], "aver") == 0)
                aver = true;
            else if (strcmp(argv[i + 1], "single") == 0)
                single = true;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--path") == 0)
        {
            mypath = argv[i + 1];
        }
        if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--start") == 0)
        {
            start_frame = atoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-lmp") == 0)
            lmptrj = true;
        if (strcmp(argv[i], "-xyz") == 0)
            xyztrj = true;
        if (strcmp(argv[i], "-mc") == 0)
            more_cols_flag = true;
    }
    if (rank == 0 && single == aver)
    {
        ierr = 2;
        cerr << "Both modes are turned on or off, please, check mode selection. Only single mode should be selected. For help use key -h or --help." << endl;
        MPI_Finalize();
        return 2;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    float **x, **y, **z;
    int bead_begin, bead_end;

    if (lmptrj)
    {
        input_file.open(mypath, ios::in);
        if (!input_file.is_open())
        {
            cerr << "Somthing went wrong with opening lmp file" << endl;
            MPI_Finalize();
            return 15;
        }
        numberFrames = 0;
        while (getline(input_file, line))
        {
            if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
                ++numberFrames;
        }
        if (start_frame >= numberFrames)
        {
            cerr << "Start frame is bigger than number of frames in file" << endl;
            MPI_Finalize();
            return 16;
        }
        else
        {
            numberFrames -= start_frame;
        }
        input_file.close();
        if (rank == 0)
            cout << "1st reading is over by " << rank << endl;
        input_file.open(mypath, ios::in);
        getline(input_file, emp);
        while (!(emp.find("ITEM: NUMBER OF ATOMS") != std::string::npos))
        {
            getline(input_file, emp);
        }
        input_file >> numberBeads;
        if (rank == 0)
            cout << "number of frames " << numberFrames << '\n'
                 << "start frame " << start_frame << '\n'
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
        while (getline(input_file, emp))
        {
            frame++;
            while (!(emp.find("ITEM: ATOMS id type xu yu zu") != std::string::npos))
            {
                if (!getline(input_file, emp))
                    break;
            }
            if (frame < start_frame)
            {
                for (int i = 0; i < numberBeads; i++)
                {
                    getline(input_file, emp);
                }
                continue;
            }
            else
            {
                for (int i = 0; i < numberBeads; i++)
                {
                    if (input_file.peek() == EOF)
                        break;
                    input_file >> id;
                    id--;
                    if (id >= bead_begin && id < bead_end)
                    {
                        if (more_cols_flag)
                        {
                            input_file >> emp >> x[id][frame - start_frame] >> y[id][frame - start_frame] >> z[id][frame - start_frame] >> emp >> emp >> emp;
                        }
                        else
                            input_file >> emp >> x[id][frame - start_frame] >> y[id][frame - start_frame] >> z[id][frame - start_frame];
                    }
                    else
                    {
                        //input_file >> emp >> emp >> emp >> emp;
                        getline(input_file, emp);
                    }
                }
            }
        }
        input_file.close();
        if (rank == 0)
            cout << "2nd reading is over by " << rank << endl;
    }
    else if (xyztrj)
    {
        input_file.open(mypath, ios::in);
        if (!input_file.is_open())
        {
            cerr << "Somthing went wrong with opening lmp file" << endl;
            MPI_Finalize();
            return 15;
        }
        numberFrames = 0;
        while (getline(input_file, line))
        {
            ++numberFrames;
        }
        numberFrames = numberFrames / 1001;
        input_file.close();
        if (rank == 0)
            cout << "1st reading is over by " << rank << endl;
        input_file.open(mypath, ios::in);
        numberBeads = 1000;
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
        while (getline(input_file, emp))
        {
            frame++;
            for (int i = 0; i < numberBeads; i++)
            {
                if (input_file.peek() == EOF)
                    break;
                if (i >= bead_begin && i < bead_end)
                {
                    input_file >> x[i][frame] >> y[i][frame] >> z[i][frame];
                }
                else
                {
                    input_file >> emp >> emp >> emp;
                }
            }
        }
        input_file.close();
        if (rank == 0)
            cout << "2nd reading is over by " << rank << endl;
    }

    if (rank == 0)
        cout << numberBeads << ' ' << reg_size << ' ' << unreg_size << ' ' << commsize << endl;
    if (aver)
    {
        float *msd = new float[numberFrames];
        float *msd_count = new float[numberFrames];
        for (int i = 0; i < numberFrames; i++)
        {
            msd[i] = 0;
            msd_count[i] = 0;
        }
        msd_aver(x, y, z, msd, msd_count, numberFrames, bead_begin, bead_end);
        float *msd_total = new float[numberFrames];
        float *msd_total_count = new float[numberFrames];
        for (int i = 0; i < numberFrames; i++)
        {
            msd_total[i] = 0;
            msd_total_count[i] = 0;
        }
        MPI_Reduce(msd, msd_total, numberFrames, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(msd_count, msd_total_count, numberFrames, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            fstream fout;
            fout.open(mypath + "_msd_a.dat", ios::out);
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
    else if (single)
    {
        int curr_bd, rcv_rnk;
        float **msd = new float *[numberBeads];
        float **msd_count = new float *[numberBeads];
        float **msd_total = new float *[numberBeads];
        float *rcv_msd = new float[numberBeads];
        float **msd_total_final = new float *[numberBeads];
        for (unsigned int i = 0; i < numberBeads; i++)
        {
            msd[i] = new float[numberFrames];
            msd_count[i] = new float[numberFrames];
            msd_total[i] = new float[numberFrames];
            msd_total_final[i] = new float[numberFrames];
            rcv_msd[i] = 0;
            for (unsigned int j = 0; j < numberFrames; j++)
            {
                msd[i][j] = 0;
                msd_count[i][j] = 0;
                msd_total[i][j] = 0;
                msd_total_final[i][j] = 0;
            }
        }
        msd_single(x, y, z, msd, msd_count, numberFrames, bead_begin, bead_end);
        for (unsigned int i = bead_begin; i < bead_end; i++)
        {
            for (unsigned int j = 0; j < numberFrames; j++)
            {
                if (msd_count[i][j] != 0)
                    msd_total[i][j] = msd[i][j] / msd_count[i][j];
                else
                    msd_total[i][j] = 0;
            }
        }

        for (unsigned int i = 0; i < numberBeads; i++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0 && (i >= bead_begin && i < bead_end))
                continue;
            if (rank == 0 || (i >= bead_begin && i < bead_end))
            {
                if (rank != 0)
                    MPI_Send(msd_total[i], numberFrames, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
                else
                {
                    rcv_rnk = (i - unreg_size) / reg_size + 1;
                    MPI_Recv(msd_total[i], numberFrames, MPI_FLOAT, rcv_rnk, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        if (rank == 0)
        {
            fstream fout;
            cout << "PATH: " << mypath + "_msd_s.dat" << endl;
            fout.open(mypath + "_msd_s.dat", ios::out);
            for (unsigned int j = 0; j < numberFrames; j++)
            {
                fout << j << '\t';
                for (unsigned int i = 0; i < numberBeads; i++)
                    fout << msd_total[i][j] << '\t';
                fout << '\n';
            }
            fout.close();
        }
    }

    cout << rank << " is done!" << endl;
    MPI_Finalize();
    return 0;
}
