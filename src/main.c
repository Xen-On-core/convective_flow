#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "utils/array.h"
#include "structures.h"
#include "calculations.h"

float epsilon = 1e-6;
#define NUMCHECK(number) \
            if (abs(number - 0) < epsilon) { \
                printf("Wrong value of %s: %f\n", #number, number); \
                exit(1); \
            }

void helper(const char *progname) {
    printf("%s starts a calculations convective flow of fluid in 2D case.\n\n", progname);
    printf("Usage:\n");
    printf("    %s [OPTION]...      \n", progname);
    printf("\nOptions:\n");
    printf("    -G, --grashof       Reynolds number\n");
    printf("    -P, --prandtl       Prandtl number\n");
    printf("    -R, --reynolds      Grashof number\n");
    printf("    -T, --time-points   number of nodes by X\n");
    printf("    -t, --t1            \n");
    printf("    -x, --x1            \n");
    printf("    -X, --x-points      number of nodes by Y\n");
    printf("    -y, --y1            \n");
    printf("    -Y, --y-points      number of nodes by time\n");
    printf("        --x-scale       \n");
    printf("        --y-scale       \n");
    printf("    -?, --help          \n");
}

int main(int argc, char *argv[]) {

    if (argc < 1)
    {
        printf("No arguments are given!");
        return 0;
    }

    static struct option long_options[] = {
        {"grashof", required_argument, NULL, 'G'},
        {"help", no_argument, NULL, '?'},
        {"prandtl", required_argument, NULL, 'P'},
        {"reynolds", required_argument, NULL, 'R'},
        {"t_1", required_argument, NULL, 't'},
        {"time-points", required_argument, NULL, 'T'},
        {"x_1", required_argument, NULL, 'x'},
        {"x-points", required_argument, NULL, 'X'},
        {"x-scale", required_argument, NULL, 1},
        {"y_1", required_argument, NULL, 'y'},
        {"y-points", required_argument, NULL, 'Y'},
        {"y-scale", required_argument, NULL, 2},
        {NULL, 0, NULL, 0}
	};
    int c;
    int option_index = 0;
    StarterPack STPKG = {};
    STPKG.Re = 1.0;
    STPKG.Gr = 10000.0;
    STPKG.Pr = 1.0;
    STPKG.N = STPKG.M = 10;
    STPKG.K = 1000;
    STPKG.y_1 = 1.0; 
    STPKG.y_1 = 1.0;
    STPKG.t_1 = 1.0;
    STPKG.T_0 = 1.0/2.0;
    STPKG.T_1 = 1.0;

    float x_scale,y_scale;

    while ((c = getopt_long(argc, argv, "G:P:R:t:T:x:X:y:Y:",
                 long_options, &option_index)) != -1)
    {
        switch (c)
        {
            case 1:
                
                x_scale = atof(optarg);
                NUMCHECK(x_scale)
                break;
            case 2:
                y_scale = atof(optarg);
                NUMCHECK(y_scale)
                break;
            case 'X':
                STPKG.N = atoi(optarg);
                if (STPKG.N < 10)
                {
                    printf("Wrong number of X steps: %d\n", STPKG.N);
                    exit(1);
                }
                break;
            case 'Y':
                STPKG.M = atoi(optarg);
                if (STPKG.M < 10)
                {
                    printf("Wrong number of Y steps: %d\n", STPKG.M);
                    exit(1);
                }
                break;
            case 'T':
                STPKG.K = atoi(optarg);
                if (STPKG.K < 10)
                {
                    printf("Wrong number of time steps: %d\n", STPKG.K);
                    exit(1);
                }
                break;
            case 'x':
                STPKG.x_1 = atof(optarg);
                NUMCHECK(STPKG.x_1)
                break;
            case 'y':
                STPKG.y_1 = atof(optarg);
                NUMCHECK(STPKG.y_1)
                break;
            case 't':
                STPKG.t_1 = atof(optarg);
                NUMCHECK(STPKG.t_1)
                break;
            case 'R':
                STPKG.Re = atof(optarg);
                NUMCHECK(STPKG.Re)
                break;
            case 'G':
                STPKG.Gr = atof(optarg);
                NUMCHECK(STPKG.Gr)
                break;
            case 'P':
                STPKG.Pr = atof(optarg);
                NUMCHECK(STPKG.Pr)
                break;
            case '?':
                helper(argv[0]);
                exit(0);
                break;
            default:
                printf("Try \"./convective_flow --help\" for more information.\n");
                exit(1);
        }
    }

   FindSolution(STPKG);

    return 0;
}
