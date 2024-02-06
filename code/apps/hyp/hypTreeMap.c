//
//  hypTreeMap.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2021.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2021.
//
//  Copyright 2019 - 2021 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the GNU Public Licence v3.0
//

#include <stdlib.h>
#include <sys/timeb.h>

#include "stopWatch.h"
#include "hypTreeMap.h"
#include "hypRaw.h"
#include "treeMap.h"
#include "hypRawCount.h"
#include "io.h"

static bool writeAndCheck(tmap utr, int per, int st, int en, bool real) {
    char mem[30], fileName[150];

    int lev = utr->maxLevel + 2;
    if(st < 0 || en < 0) {
        snprintf(fileName, 149, real ? "%s/%02d/hypReal%02d_res%02d.tmap" :
                 "%s/%02d/hyp%02d_res%02d.tmap", HYP_FOLDER, per, per, lev);
    } else {
        snprintf(fileName, 149, real ? "%s/%02d/hypReal%02d_res%02d_%d-%d.tmap" :
                 "%s/%02d/hyp%02d_res%02d_%d-%d.tmap", HYP_FOLDER, per, per, lev, st, en);
    }
    
    mem_size(tmap_memory(utr), mem);
    if(real) {
        printf("\nThe real tree of resolution %d uses %s of RAM to store %ld points.\n",
               utr->maxLevel + 2, mem, utr->count);
    } else {
        printf("\nThe tree of resolution %d uses %s of RAM to store %ld points.\n",
               utr->maxLevel + 2, mem, utr->count);
    }

    if(utr->count == 0) {
        printf("Empty maps are not stored.\n");
        
        return true;
    }
    
    bool ok = tmap_save(utr, fileName);
    if(ok) {
        file_size(io_size_of(fileName), mem);
        
        printf("Points map written to %s, of size %s.\n", fileName, mem);
        
        tmap ctr = tmap_load(fileName, 0, true, NULL);
        if(ctr != NULL && tmap_eq(utr, ctr)) {
            printf("The file is ready for publishing.\n");
        } else {
            printf("The file is CORRUPTED !!!\n");
            
            ok = false;
        }
        
        if(ctr != NULL) {
            tmap_free(ctr);
        }
    } else {
        printf("Could not write the treeMap to %s.\n", fileName);
        
        ok = false;
    }
    
    fflush(stdout);
    
    return ok;
}

tmap mapReal(nset ps, int mlev, int lstep) {
    if(ps == NULL) {
        return NULL;
    }
    
    tmap tr = tmap_new(0, 0, -2, mlev, lstep, -2);
    
    u128 p, p0;
    p0->y = 0;
    ulong eps = ps->eps;
    for (ulong i = 0; i < ps->count; i++) {
        nset_point(p, ps, i);
        p0->x = p->x;
        
        if(u128_eq(p, p0, eps)) {
            tmap_add(tr, p0);
        }
    }
    
    return tr;
}

static bool save(tmap tree, tmap treal, int per, int st, int en, int minLev, int levStep) {
    bool ok = writeAndCheck(tree, per, st, en, false);
    ok = ok && writeAndCheck(treal, per, st, en, true);
    
    while(ok && minLev < tree->maxLevel + levStep) {
        ok = tmap_reduce(tree, tree->maxLevel - levStep);
        ok = ok && tmap_reduce(treal, treal->maxLevel - levStep);
        
        if(ok) {
            ok = ok && writeAndCheck(tree, per, st, en, false);
            ok = ok && writeAndCheck(treal, per, st, en, true);
        } else {
            printf("Could not reduce the resolution of the map.\n");
        }
    }
    
    if(! ok) {
        printf("\n    ****\nSome error occurred, please check the messages above !\n    ****\n");
    } else {
        printf("\n");
    }
    
    fflush(stdout);

    return ok;
}

bool map(int per, int totLev, int minLev, int maxLev, int levStep, int group, int start, int end) {
    if(totLev > minLev || totLev < 10 || minLev > maxLev || maxLev > TMAP_MAX_LEVEL ||
       levStep < 2 || levStep > TMAP_MAX_LEVEL_STEP || (maxLev - minLev) % levStep != 0 ||
       (maxLev - totLev) % levStep != 0 || group < 1) {
        printf("The parameters to map() are not valid, please try again.\n");
        
        return false;
    }
    
    if(minLev == maxLev) {
        printf("Creating maps of hyperbolic centers of period %d, with resolution per unit 2^%d\n",
               per, maxLev);
    } else if(start < 0 || end < 0) {
        printf("Creating maps of hyperbolic centers of period %d, with resolutions per unit from 2^%d to 2^%d\n",
               per, minLev, maxLev);
    } else {
        printf("Creating maps of hyperbolic centers of period %d, with resolutions per unit from 2^%d to 2^%d,\nfrom nset files with indexes from %d to %d.\n",
               per, minLev, maxLev, start, end - 1);
    }
    
    fflush(stdout);

    char time[80], t1[80], t2[80];
    struct timeb ats, ts;
    ftime(&ats);
    
    time_stamp(time, 80, true, true);
    printf("Started at %s\n\n", time);
    fflush(stdout);

    char fileName[100];
    bool ok = true;
    tmap utr, urt, ctr, rtr, ttr, trt;
    utr = tmap_new_mandel(maxLev - levStep, levStep);
    urt = tmap_new_mandel(maxLev - levStep, levStep);
    ttr = tmap_new_mandel(totLev - levStep, levStep);
    trt = tmap_new_mandel(totLev - levStep, levStep);
    
    // check the range
    int files = hyp_resultsCount(per);
    if(start < 0 || end < 0 || start > end) {
        start = 0;
        end = files;
    } else {
        end = end > files ? files : end;
    }
    
    // load all final nset files and create their maps
    for (int i = start; ok && i < end; i++) {
        ftime(&ts);
        
        nset ps = hyp_loadResults(per, i);
        ctr = ps == NULL ? NULL : tmap_map(ps, maxLev - levStep, levStep);
        rtr = mapReal(ps, maxLev - levStep, levStep);
        nset_free(ps);

        hyp_fileName(fileName, 99, per, i);
        if(ctr == NULL || rtr == NULL) {
            printf("Failed to load %s\n", fileName);
            fflush(stdout);
            
            ok = false;
        } else {
            lapse(&ts, time);

            ftime(&ts);
            tmap_union(utr, ctr);
            tmap_union(urt, rtr);
            lapse(&ts, t1);
            
            ftime(&ts);
            tmap_free(ctr);
            tmap_free(rtr);
            lapse(&ts, t2);
            
            printf("Loaded and mapped %s in %s (union in %s, free in %s)\n", fileName, time, t1, t2);
            fflush(stdout);
            
            // save every group nset files
            if((i - start + 1) % group == 0) {
                ftime(&ts);
                
                int st = i - group + 1;
                st = st == 0 && i == files - 1 ? -1 : st; // check if all files are in this map
                
                ok = ok && save(utr, urt, per, st, i, minLev, levStep);
                
                ok = ok && tmap_reduce(utr, totLev - levStep);
                ok = ok && tmap_union(ttr, utr);
                tmap_free(utr);
                
                ok = ok && tmap_reduce(urt, totLev - levStep);
                ok = ok && tmap_union(trt, urt);
                tmap_free(urt);

                utr = tmap_new_mandel(maxLev - levStep, levStep);
                urt = tmap_new_mandel(maxLev - levStep, levStep);

                lapse(&ts, time);
                printf("Maps reduced, saved, read, compared and freed in %s\n\n", time);
                fflush(stdout);
            }
        }
    }
    
    // if the number of files is not a multiple of @c group
    if(ok && utr->count > 0) {
        ftime(&ts);
        
        int f = end - 1;
        int prev = group * ((f - start) / group) + start;
        prev = prev == 0 && end == files ? -1 : prev; // check if all files are in this map
        ok = ok && save(utr, urt, per, prev, end - 1, minLev, levStep);
        
        ok = ok && tmap_reduce(utr, totLev - levStep);
        ok = ok && tmap_union(ttr, utr);
        tmap_free(utr);
        
        ok = ok && tmap_reduce(urt, totLev - levStep);
        ok = ok && tmap_union(trt, urt);
        tmap_free(urt);

        lapse(&ts, time);
        printf("Maps reduced and saved in %s\n", time);
        fflush(stdout);
    } else {
        tmap_free(utr);
        tmap_free(urt);
    }
    
    // save the global maps of level totLev
    ftime(&ts);
    int s = start == 0 && end == files ? 0 : start;
    save(ttr, trt, per, s, end - 1, minLev, levStep);
    lapse(&ts, time);
    
    printf("Global maps saved in %s\n", time);
    fflush(stdout);
    
    tmap_free(ttr);
    tmap_free(trt);
    
    lapse(&ats, time);
    printf("Total running time %s\n", time);
    
    return ok;
}

// MARK: the help system and the main function

static const char* before = "This task maps the hyperbolic centers of given period.\n\n";
static const char* after = "\nSeveral resolutions of maps can be produced and stored in one pass.\n\n";

static const char *parameters[] = {
    "period",
    "group",
    "first file",
    "last file + 1"
};

static const char *types[] = {
    "required",
    "optional",
    "optional",
    "optional"
};

static const char *defaults[] = {
    "",
    "8",
    "0",
    "file count"
};

static const char *descriptions[] = {
    "the period, integer, at least 3, at most 41",
    "the number of nset files in a map of resolution 2^30",
    "the index of the first nset file to use",
    "the index of the first nset file to ignore"
};

static const char *headers[] = {
    "Parameter",
    "Type",
    "Default value",
    "Description"
};

static const int paramCount = 4;
static const int columnWidths[] = {18, 18, 18};

/// Prints instructions for usage and some details about the command line arguments.
static void help(void) {
    printf("%s", before);
    char format[50];
    snprintf(format, 45, "    %%-%ds %%-%ds %%-%ds %%s\n",
             columnWidths[0], columnWidths[1], columnWidths[2]);
    
    printf(format, headers[0], headers[1], headers[2], headers[3]);
    printf("\n");
    for(int i = 0; i < paramCount; i++) {
        printf(format, parameters[i], types[i], defaults[i], descriptions[i]);
    }
    
    printf("%s", after);
    fflush(stdout);
}

int hyp_tree_main(int argc, const char * argv[]) {
    int per, gr = HYP_TREE_GROUP, st = -1, en = -1;
    
    if(argc < 1 || sscanf(argv[0], "%d", &per) < 1 || per < 1 || per > HYP_RAW_MAX_PER) {
        help();
        
        return 1;
    }
    
    if(argc >= 2 && (sscanf(argv[1], "%d", &gr) < 1)) {
        help();
        
        return 1;
    }
    
    if(argc >= 4 && (sscanf(argv[2], "%d", &st) < 1 || sscanf(argv[3], "%d", &en) < 1 ||
       st < 0 || en <= st)) {
        help();
        
        return 1;
    }
    
    return ! map(per, HYP_TREE_GLOBAL_LEV, HYP_TREE_MIN_LEV, HYP_TREE_MAX_LEV, HYP_TREE_LEV_STEP, gr, st, en);
}
