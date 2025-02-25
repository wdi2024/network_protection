#include "proc_trav.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = wpalloc(1,sizeof(PreprocData));
    //no mode 2, so it is just randomly picking up children
    pd->P = 100; // P% percent of total passes in mode 1, the remaining in mode 2, 300 is mode 3
    pd->total = 10; //number of total passes
    pd->K = 5;
    // pd->LEVEL_NUM = 50; 


    inputGraph(pd); //load graph data from standard input

    double tm = timer();
    initializeAllData(pd); //init data structure
    double tm1 = timer() - tm;
    tm = timer();

    preProcessAll(pd); // preproc by traversing the graph for pd->total times
    double tm2 = timer() - tm;
    tm = timer();

    //transFormAndSingleCut(pd); //extra weak graph, find all possible Single Cuts, and level the only vomponent

    spanningTreeBuildAll(pd);

    double tm3 = timer() - tm;

    printf("the final total time: %lf (%lf %lf %lf)\n",(tm1+tm2+tm3),tm1,tm2,tm3);
    exit(0);

    
}
