//
// Created by facundo on 11/6/20.
//

#ifndef COVID_19_STATES_H
#define COVID_19_STATES_H

#endif //COVID_19_STATES_H

enum State {freeCell, susceptible, sickNoContagion, sickContagion, isolatedSick, dead, cured};
enum Age {kid, adult, elder};
enum Sex {woman, man};

typedef struct Person{
    enum State state;
    enum Age age;
    int riskDisease;
    int riskProfession;
    int vaccines;
    int days;
    enum Sex sex;
} Person;
