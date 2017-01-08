//
//  util.c
//  SingleChain
//
//  Created by Linling Miao on 12/18/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* this file only used to parse input arguments */

#include <getopt.h>
#include "util.h"

void ParseInput(int argc, char * argv[]){
    /* Default Values */
    flowType = 0;
    flowRate = 0;
    numThread = 1;
    N = 20;
    
    int option = 0;
    
    if(argc < 9) printf("Will use default values if an input is not specified\n");
    while((option = getopt(argc, argv, "n:t:f:r:")) != -1){
        switch(option){
            case 'n':
                sscanf(optarg, "%d", &N);
                break;
            case 't':
                sscanf(optarg, "%d", &numThread);
                break;
            case 'f':
                sscanf(optarg, "%d", &flowType);
                break;
            case 'r':
                flowratec = optarg;
                sscanf(optarg, "%lf", &flowRate);
                break;
            case '?':
                printf("Unknown option -%c.\n Execution abort...", optopt);
                exit(EXIT_FAILURE);
        }
    }
    
    printf("Chain length: %d\n", N);
    printf("Number of chains: %d\n", NP);
    printf("Box size: %f\n",L);
    printf("Flow type: %d, Flow rate: %lf\n", flowType,flowRate);
    printf("Running with number of %d OpenMP threads...\n", numThread);
}
