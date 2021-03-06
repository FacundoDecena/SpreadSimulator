#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "person.h"

#define FALSE 0
#define TRUE 1
#define DAYS 120                // numbers of days
#define DENSITY 0.50            // population density
#define INITIAL_INFECTION 0.002 // % of initial cases
#define KIDS 0.30               // % of the population
#define ADULTS 0.54             // % of the population
//      ELDERS 0.16             // % of the population

#define ILLNESS_STRENGTH 2.4     // enhancers
#define KID_SUSCEPTIBILITY 0.3   // enhancers
#define ADULT_SUSCEPTIBILITY 0.5 // enhancers
#define ELDER_SUSCEPTIBILITY 0.9 // enhancers
#define RISK_FACTOR 0.15         // enhancers
#define VACCINES_BOOST 0.005     // enhancers
#define SIZE 1500                // Matrix size for rows and columns

// Returns a initialize person
Person generatePerson();

// Adds people to the matrix
void init(Person *matrix);

// susToSick changes the state of a susceptible person to sick without contagion if
// meets the conditions and returns TRUE if it changed, FALSE if it did not.
// *matrix is a vector[9] (representing a 3x3 matrix) where de susceptible person is in the position 4
int susToSick(Person *matrix, Person *p);

// noConToCon changes the state of a sick person without contagion to a contagious one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
int noConToCon(Person *p);

// conToAis changes the state of a sick person with contagion to a isolated one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
int conToAis(Person *p);

// conToAis changes the state of a sick person with contagion to a isolated one if
// meets the conditions. returns TRUE if it changed, FALSE if it did not.
int sickToD_or_R(Person *p);

// returns an array with de indexes of the neighbors of a given index
void neighbors(int index, int *n);

int main()
{
    int i;
    Person *matrix = (Person *)malloc(sizeof(Person) * SIZE * SIZE);
    Person *lastMatrix = (Person *)malloc(sizeof(Person) * SIZE * SIZE);
    Person *neighborhood = (Person *)malloc(sizeof(Person) * 9);

    int indexes[9];
    init(lastMatrix);
   
    srand(time(0));
    double start, end;
    double cpu_time_used, sum = 0.0;
    /*
    for (int k = 0; k < SIZE * SIZE; ++k) {
        printf("%d ", lastMatrix[k].state);
                if ((k+1) % SIZE == 0)
                    printf("\n");
    }
    printf("\n");
    */
    for (int numRep = 0; numRep < 10; numRep++)
    {

        start = omp_get_wtime();
        for (i = 0; i < DAYS; i++)
        {
            
    #pragma omp parallel for schedule(static, 6)
            for (int j = 0; j < SIZE * SIZE; j++)
            {
                if (lastMatrix[j].state == freeCell)
                    continue;
                neighbors(j, indexes);
                neighborhood[0] = lastMatrix[indexes[0]];
                neighborhood[1] = lastMatrix[indexes[1]];
                neighborhood[2] = lastMatrix[indexes[2]];
                neighborhood[3] = lastMatrix[indexes[3]];
                neighborhood[4] = lastMatrix[j];
                neighborhood[5] = lastMatrix[indexes[5]];
                neighborhood[6] = lastMatrix[indexes[6]];
                neighborhood[7] = lastMatrix[indexes[7]];
                neighborhood[8] = lastMatrix[indexes[8]];
                switch (neighborhood[4].state)
                {
                case susceptible:
                    susToSick(neighborhood, &neighborhood[4]);
                    break;
                case sickNoContagion:
                    neighborhood[4].days++;
                    noConToCon(&neighborhood[4]);
                    break;
                case sickContagion:
                    neighborhood[4].days++;
                    conToAis(&neighborhood[4]);
                    break;
                case isolatedSick:
                    neighborhood[4].days++;
                    break;
                default:
                    break;
                }
                sickToD_or_R(&neighborhood[4]);
                matrix[j] = neighborhood[4];
            }
            memcpy(lastMatrix, matrix, (SIZE * SIZE) * sizeof(Person));
        }
        end = omp_get_wtime();
        cpu_time_used = (end - start);
        sum += cpu_time_used;
    }
    /*
    printf("\n");
    for (int k = 0; k < SIZE * SIZE; ++k) {
        printf("%d ", lastMatrix[k].state);
        if ((k+1) % SIZE == 0)
            printf("\n");
    }
    */
    printf("Tiempo promedio: %lf\n\n", sum/10);
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

    ran = (double)rand() / (double)RAND_MAX;
    if (ran < KIDS)
        p.age = kid;
    else if (ran < KIDS + ADULTS)
        p.age = adult;
    else
        p.age = elder;

    ran = (double)rand() / (double)RAND_MAX;
    if (ran < 0.5)
        p.riskDisease = FALSE;
    else
        p.riskDisease = TRUE;

    ran = (double)rand() / (double)RAND_MAX;
    if (ran < 0.5)
        p.riskProfession = FALSE;
    else
        p.riskProfession = TRUE;

    ran = (double)rand() / (double)RAND_MAX;
    if (ran < 0.5)
        p.sex = man;
    else
        p.sex = woman;

    ran = (double)rand() / (double)RAND_MAX;
    if (ran < 0.5)
        p.vaccines = TRUE;
    else
        p.vaccines = FALSE;

    ran = (double)rand() / (double)RAND_MAX;
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

void neighbors(int index, int *n)
{
    n[0] = index - SIZE - 1;
    n[1] = index - SIZE;
    n[2] = index - SIZE + 1;
    n[3] = index - 1;
    n[5] = index + 1;
    n[6] = index + SIZE - 1;
    n[7] = index + SIZE;
    n[8] = index + SIZE + 1;
    //Lados
    if (index % SIZE == 0)
    {
        if (index == 0)
        {
            n[0] = SIZE * SIZE - 1;
            n[1] = SIZE * SIZE - SIZE;
            n[2] = SIZE * SIZE - SIZE + 1;
            n[3] = SIZE - 1;
            n[6] = 2 * SIZE - 1;
            return;
        }
        if (index == SIZE * SIZE - SIZE)
        {
            n[0] = index - 1;
            n[3] = SIZE * SIZE - 1;
            n[6] = SIZE - 1;
            n[7] = 0;
            n[8] = 1;
            return;
        }
        n[0] = index - 1;
        n[3] = index + SIZE - 1;
        n[6] = index + 2 * SIZE - 1;
        return;
    }
    if ((index + 1) % SIZE == 0)
    {
        if (index == SIZE - 1)
        {
            n[0] = SIZE * SIZE - 2;
            n[1] = SIZE * SIZE - 1;
            n[2] = SIZE * SIZE - SIZE;
            n[5] = 0;
            n[8] = index + 1;
            return;
        }
        if (index == SIZE * SIZE - 1)
        {
            n[2] = index - 2 * SIZE + 1;
            n[5] = SIZE * SIZE - SIZE;
            n[6] = SIZE - 2;
            n[7] = SIZE - 1;
            n[8] = 0;
            return;
        }
        n[2] = index + 1 - 2 * SIZE;
        n[5] = index + 1 - SIZE;
        n[8] = index + 1;
        return;
    }
    if (0 < index && index < SIZE - 1)
    {
        n[0] = index - 1 + SIZE * SIZE - SIZE;
        n[1] = index + 0 + SIZE * SIZE - SIZE;
        n[2] = index + 1 + SIZE * SIZE - SIZE;
        return;
    }
    if (SIZE * SIZE - SIZE < index && index < SIZE * SIZE - 1)
    {
        n[6] = index - 1 - SIZE * SIZE + SIZE;
        n[7] = index + 0 - SIZE * SIZE + SIZE;
        n[8] = index + 1 - SIZE * SIZE + SIZE;
        return;
    }
}

int susToSick(Person *matrix, Person *p)
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
        if (matrix[i].state == sickContagion)
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

    if ((double)rand() / (double)RAND_MAX < probContagion)
    {
        p->state = sickNoContagion;
        return TRUE;
    }
    return FALSE;
}

int noConToCon(Person *p)
{
    if (p->days == 4)
    {
        p->state = sickContagion;
        return TRUE;
    }
    return FALSE;
}

int conToAis(Person *p)
{
    double ran = (double)rand() / (double)RAND_MAX;
    if (p->days == 7)
    {
        if (ran > 0.1)
        {
            p->state = isolatedSick;
            return TRUE;
        }
    }
    return FALSE;
}

int sickToD_or_R(Person *p)
{
    double ran = (double)rand() / (double)RAND_MAX;
    if (p->days == 14)
    {
        p->days = 0;
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
        if (ran < probDead)
        {
            p->state = dead;
            return TRUE;
        }
        else
        {
            p->state = cured;
            return FALSE;
        }
    }
    return -1;
}
