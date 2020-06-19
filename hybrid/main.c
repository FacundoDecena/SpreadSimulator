#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include "person.h"

#define MASTER 0
#define FALSE 0
#define TRUE 1
#define DAYS 120                // numbers of days
#define DENSITY 0.50            // population density
#define INITIAL_INFECTION 0.002 // % of initial cases
#define KIDS 0.30               // % of the population
#define ADULTS 0.54             // % of the population
//      ELDERS 0.16                 // % of the population

#define ILLNESS_STRENGTH 2.4     // enhancers
#define KID_SUSCEPTIBILITY 0.3   // enhancers
#define ADULT_SUSCEPTIBILITY 0.5 // enhancers
#define ELDER_SUSCEPTIBILITY 0.9 // enhancers
#define RISK_FACTOR 0.15         // enhancers
#define VACCINES_BOOST 0.005     // enhancers
#define SIZE 24                  // Matrix size for rows and columns

// Returns a initialize person
Person generatePerson();

// Adds people to the matrix
void init(Person *matrix);

// susToSick changes the state of a susceptible person to sick without contagion if
// meets the conditions and returns TRUE if it changed, FALSE if it did not.
// *matrix is a vector[9] (representing a 3x3 matrix) where de susceptible person is in the position 4
int susToSick(Person **matrix, Person *p);

// noConToCon changes the state of a sick person without contagion to a contagious one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
void noConToCon(Person *p);

// conToAis changes the state of a sick person with contagion to a isolated one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
void conToAis(Person *p);

// conToAis changes the state of a sick person with contagion to a isolated one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
void sickToD_or_R(Person *p);

// returns an array with de indexes of the neighbors of a given index
//void neighbors(int index, int *n);
void neighbors(Person *matrix, int matrix_h, int cell_x, int cell_y, Person **out_buffer);

int main()
{
    int i, j, k;
    double start, end;
    double cpu_time_used;
    int nprocs, myRank;

    srandom(time(0));

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (myRank == MASTER)
    {
        if (SIZE % nprocs != 0)
        {
            fprintf(stderr, "[ERR] rows (%d) %% nprocs (%d) != 0 (%d)\n", SIZE, nprocs, SIZE % nprocs);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    int rowsPerProcessor = (SIZE / nprocs) + 2;
    Person *myMatrix = (Person *)malloc(sizeof(Person) * rowsPerProcessor * SIZE);
    Person *myLastMatrix = (Person *)malloc(sizeof(Person) * rowsPerProcessor * SIZE);
    Person *matrix;
    Person *lastMatrix;
    int *sendCounts;
    int *recvCounts;
    int *displacements;

    MPI_Datatype MPI_PERSON;
    int mpiCellBlockLengths[] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Aint mpiCellDisplacements[] = {
        offsetof(Person, state),
        offsetof(Person, age),
        offsetof(Person, riskDisease),
        offsetof(Person, riskProfession),
        offsetof(Person, vaccines),
        offsetof(Person, days),
        offsetof(Person, sex)};
    MPI_Datatype mpiCellLengths[] = {
        MPI_INT, // state
        MPI_INT, // age
        MPI_INT, // riskDisease
        MPI_INT, // riskProfession
        MPI_INT, // vaccines
        MPI_INT, // days
        MPI_INT  // sex
    };
    MPI_Type_create_struct(
        7,
        mpiCellBlockLengths,
        mpiCellDisplacements,
        mpiCellLengths,
        &MPI_PERSON);
    MPI_Type_commit(&MPI_PERSON);
    if (myRank == MASTER)
    {
        matrix = (Person *)malloc((size_t)(SIZE + 2) * SIZE * sizeof(Person));
        lastMatrix = (Person *)malloc((size_t)(SIZE + 2) * SIZE * sizeof(Person));
        init(&lastMatrix[SIZE]);
        memcpy(lastMatrix, &lastMatrix[SIZE * SIZE], (size_t)SIZE * sizeof(Person));
        memcpy(&lastMatrix[(SIZE + 1) * SIZE], &lastMatrix[SIZE], (size_t)SIZE * sizeof(Person));
        sendCounts = malloc((size_t)nprocs * sizeof(int));
        recvCounts = malloc((size_t)nprocs * sizeof(int));
        displacements = malloc((size_t)nprocs * sizeof(int));
        for (int i = 0; i < nprocs; i++)
        {
            sendCounts[i] = rowsPerProcessor * SIZE;
            recvCounts[i] = (rowsPerProcessor - 2) * SIZE;
            displacements[i] = i * (rowsPerProcessor - 2) * SIZE;
        }
    }
    Person *myNeighborhood[9];
    if (myRank == MASTER)
    {
        /*/Uncomment to print
        for (int index = 0; index < SIZE * SIZE; ++index) {
            printf("%d ", lastMatrix[index].state);
                    if ((index + 1) % SIZE == 0)
                        printf("\n");
        }
        printf("\n");
        /*/
        start = omp_get_wtime();
    }
    for (i = 0; i < DAYS; i++)
    {
        MPI_Scatterv(
            lastMatrix,
            sendCounts,
            displacements,
            MPI_PERSON,
            myLastMatrix,
            rowsPerProcessor * SIZE,
            MPI_PERSON,
            0,
            MPI_COMM_WORLD);

        memcpy(myMatrix, myLastMatrix, (size_t)(rowsPerProcessor * SIZE) * sizeof(Person));

#pragma omp parallel for collapse(2)
        for (k = 0; k < rowsPerProcessor - 1; ++k)
        {
            for (j = 0; j < SIZE; j++)
            {
                Person *current = &myMatrix[SIZE + (k * SIZE) + j];
                switch (current->state)
                {
                case freeCell:
                    continue;
                case susceptible:
                    neighbors(myLastMatrix, rowsPerProcessor, j, k + 1, myNeighborhood);
                    susToSick(myNeighborhood, current);
                    break;
                case sickNoContagion:
                    current->days++;
                    noConToCon(current);
                    break;
                case sickContagion:
                    current->days++;
                    conToAis(current);
                    break;
                case isolatedSick:
                    current->days++;
                    break;
                default:
                    break;
                }
                sickToD_or_R(current);
            }
        }

        MPI_Gatherv(
            &myMatrix[SIZE],
            (rowsPerProcessor - 2) * SIZE,
            MPI_PERSON,
            &matrix[SIZE],
            recvCounts,
            displacements,
            MPI_PERSON,
            0,
            MPI_COMM_WORLD);

        void *myTemp = myLastMatrix;
        myLastMatrix = myMatrix;
        myMatrix = myTemp;

        if (myRank == MASTER)
        {
            void *temp = lastMatrix;
            lastMatrix = matrix;
            matrix = temp;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myRank == MASTER)
    {
        end = omp_get_wtime();
        cpu_time_used = end - start;
        /*/Uncomment to print
        for (int index = 0; index < SIZE * SIZE; ++index) {
            printf("%d ", matrix[index].state);
            if ((index + 1) % SIZE == 0)
                printf("\n");
        }
        printf("\n");
        /*/
        printf("Tiempo transcurrido: %lf\n\n", cpu_time_used);

        free(matrix);
        free(lastMatrix);
        free(sendCounts);
        free(recvCounts);
        free(displacements);
    }
    free(myMatrix);
    free(myLastMatrix);

    MPI_Finalize();
    return 0;
}

Person generatePerson()
{
    Person p;
    double ran = (double)rand() / (double)RAND_MAX;
    if (ran > DENSITY)
    {
        p.state = freeCell;
        return p;
    }

    ran = (double)random() / (double)RAND_MAX;
    if (ran < KIDS)
        p.age = kid;
    else if (ran < KIDS + ADULTS)
        p.age = adult;
    else
        p.age = elder;

    ran = (double)random() / (double)RAND_MAX;
    if (ran < 0.5)
        p.riskDisease = FALSE;
    else
        p.riskDisease = TRUE;

    ran = (double)random() / (double)RAND_MAX;
    if (ran < 0.5)
        p.riskProfession = FALSE;
    else
        p.riskProfession = TRUE;

    ran = (double)random() / (double)RAND_MAX;
    if (ran < 0.5)
        p.sex = man;
    else
        p.sex = woman;

    ran = (double)random() / (double)RAND_MAX;
    if (ran < 0.5)
        p.vaccines = TRUE;
    else
        p.vaccines = FALSE;

    ran = (double)random() / (double)RAND_MAX;
    if (ran < INITIAL_INFECTION)
        p.state = sickContagion;
    else
        p.state = susceptible;
    p.days = 0;

    return p;
}

void init(Person *matrix)
{
    int i;
    for (i = 0; i < (SIZE * SIZE); i++)
    {
        Person p = generatePerson();
        matrix[i] = p;
    }
}

void neighbors(Person *matrix, int matrix_h, int cell_x, int cell_y, Person **out_buffer)
{
    out_buffer[0] = &matrix[((cell_y - 1 + matrix_h) % matrix_h) * SIZE + (cell_x - 1 + SIZE) % SIZE];
    out_buffer[1] = &matrix[((cell_y - 1 + matrix_h) % matrix_h) * SIZE + (cell_x - 0 + SIZE) % SIZE];
    out_buffer[2] = &matrix[((cell_y - 1 + matrix_h) % matrix_h) * SIZE + (cell_x + 1 + SIZE) % SIZE];
    out_buffer[3] = &matrix[((cell_y - 0 + matrix_h) % matrix_h) * SIZE + (cell_x - 1 + SIZE) % SIZE];
    out_buffer[3] = &matrix[((cell_y - 0 + matrix_h) % matrix_h) * SIZE + (cell_x - 0 + SIZE) % SIZE];
    out_buffer[5] = &matrix[((cell_y + 0 + matrix_h) % matrix_h) * SIZE + (cell_x + 1 + SIZE) % SIZE];
    out_buffer[6] = &matrix[((cell_y + 1 + matrix_h) % matrix_h) * SIZE + (cell_x - 1 + SIZE) % SIZE];
    out_buffer[7] = &matrix[((cell_y + 1 + matrix_h) % matrix_h) * SIZE + (cell_x + 0 + SIZE) % SIZE];
    out_buffer[8] = &matrix[((cell_y + 1 + matrix_h) % matrix_h) * SIZE + (cell_x + 1 + SIZE) % SIZE];
}

int susToSick(Person **matrix, Person *p)
{
    int i, sickNeighbors = 0;
    double probContagion = 0.0, ageSusceptibility = 0.0;

    switch (p->age)
    {
    case kid:
        ageSusceptibility = KID_SUSCEPTIBILITY;
        break;
    case adult:
        ageSusceptibility = ADULT_SUSCEPTIBILITY;
        break;
    case elder:
        ageSusceptibility = ELDER_SUSCEPTIBILITY;
        break;
    }

    for (i = 0; i < 9; i++)
    {
        if (i == 5)
            continue;
        if (matrix[i]->state == sickContagion)
            sickNeighbors += 1;
    }
    if (sickNeighbors != 0)
    {
        probContagion = ((double)sickNeighbors / 8 * ILLNESS_STRENGTH) + ageSusceptibility;

        if (p->riskDisease)
            probContagion += RISK_FACTOR;

        if (p->riskProfession)
            probContagion += RISK_FACTOR;

        probContagion /= 7;
    }

    if ((double)random() / (double)RAND_MAX < probContagion)
    {
        p->state = sickNoContagion;
        return TRUE;
    }
    return FALSE;
}

void noConToCon(Person *p)
{
    if (p->days == 4)
        p->state = sickContagion;
}

void conToAis(Person *p)
{
    double rand = (double)random() / (double)RAND_MAX;
    if (p->days == 7)
    {
        if (rand > 0.1)
            p->state = isolatedSick;
    }
}

void sickToD_or_R(Person *p)
{
    double rand = (double)random() / (double)RAND_MAX;
    if (p->days == 14)
    {
        double probDead = 0.0;
        switch (p->age)
        {
        case kid:
            probDead = 0.01;
            break;
        case adult:
            probDead = 0.013;
            break;
        case elder:
            probDead = 0.148;
        }
        if (p->vaccines)
            probDead -= VACCINES_BOOST;
        if (rand < probDead)
            p->state = dead;
        else
            p->state = cured;
        p->days = 0;
    }
}
