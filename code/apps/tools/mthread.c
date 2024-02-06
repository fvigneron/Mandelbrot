//
//  mthread.c
//
//  Authors: Nicolae Mihalache, François Vigneron
//
//  Please cite the references below if you use or distribute this software.
//
//  • [1] N. Mihalache & F. Vigneron. How to split a tera-polynomial. 2022.
//  • [2] N. Mihalache & F. Vigneron. How to compute the coefficients of holomorphic maps. 2022.
//
//  Copyright 2019 - 2022 Univ. Paris-Est Créteil, Univ. de Reims Champagne-Ardenne.
//  This software is released under the Lesser GNU Public Licence v3.0
//

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>

#include "stopWatch.h"
#include "mthread.h"
#include "ups.h"
#include "io.h"

static mth_threads mth_all = {0, 0, NULL};

static inline bool mth_lock(void) {
    if(pthread_mutex_lock(&mth_all.lock) != 0) {
        io_error_echo("\n#### #### Lockdead ! #### ####\n\n");
        
        return false;
    }
    
    return true;
}

static inline bool mth_unlock(void) {
    if(pthread_mutex_unlock(&mth_all.lock) != 0) {
        io_error_echo("\n#### #### Deadlock ! #### ####\n\n");
        
        return false;
    }
    
    return true;
}

static ulong get_key(mthread *task) {
    if(! mth_lock()) {
        return -1;
    }
    
    ulong key = task->key;
    
    return mth_unlock() ? key : -1;
}

static ulong get_state(mthread *task) {
    if(! mth_lock()) {
        return -1;
    }
    
    ulong state = task->state;
    
    return mth_unlock() ? state : -1;
}

static bool set_state(mthread *task, byte state) {
    if(! mth_lock()) {
        return false;
    }
    
    task->state = state;
    
    return mth_unlock();
}

static bool set_result(mthread *task, void *result) {
    if(! mth_lock()) {
        return false;
    }
    
    task->result = result;
    task->state = result == NULL ? MTH_STATE_ERROR : MTH_STATE_DONE;
    
    return mth_unlock();
}

static void *mth_std_run(void *ta) {
    mthread *task = ta;
    struct timespec sleep = {.tv_sec = task->nano / BILLION, .tv_nsec = task->nano % BILLION};
    
    while(get_key(task) != MTH_KEY); // double check
    
    while(get_state(task) != MTH_STATE_EXIT) {
        if(! mth_lock()) {
            task->result = NULL;
            task->state = MTH_STATE_ENDED;
            
            return NULL;
        }
        
        if(task->state == MTH_STATE_READY && task->run != NULL && task->task != NULL) {
            task->state = MTH_STATE_RUNNING;
            
            if(mth_unlock()) {
                set_result(task, task->run((void *) task->task, task->ID, task->index));
            } else {
                task->result = NULL;
                task->state = MTH_STATE_ENDED;
                
                return NULL;
            }
        } else {
            if(mth_unlock()) {
                nanosleep(&sleep, NULL);
            } else {
                task->result = NULL;
                task->state = MTH_STATE_ENDED;
                
                return NULL;
            }
        }
    }
    
    set_state(task, MTH_STATE_ENDED);
    
    return ta;
}

uint mth_threads_count(void) {
    return mth_all.count;
}

bool mth_init(uint threads, mth_run run, ulong us_thread, ulong us_main) {
    if(mth_all.count > 0 || run == NULL || threads == 0 || threads > MTH_MAX_THREADS) {
        return false;
    }
    
    if(pthread_mutex_init(&mth_all.lock, NULL) != 0) {
        return false;
    }
    
    mth_all.threads = malloc(sizeof(mthread) * threads);
    if(mth_all.threads == NULL) {
        return false;
    }
    
    for (int i = 0; i < threads; i++) {
        mth_all.threads[i].state = MTH_STATE_ERROR;
        mth_all.threads[i].index = -1;
        mth_all.threads[i].key = 0;
    }
    
    bool ok = true;
    for (int i = 0; i < threads && ok; i++) {
        mth_all.threads[i].state = MTH_STATE_WAIT;
        mth_all.threads[i].index = i;
        mth_all.threads[i].nano = us_thread * 1000;
        mth_all.threads[i].run = run;
        mth_all.threads[i].task = NULL;
        mth_all.threads[i].result = NULL;
        mth_all.threads[i].ID = -1;
        
        ok = ok && pthread_create((pthread_t *) &mth_all.threads[i].thread, NULL, &mth_std_run,
                                  (void *) (mth_all.threads + i)) == 0;
        
        mth_all.threads[i].key = MTH_KEY;
    }
    
    if(! ok) {
        for (int i = 0; i < threads && ok; i++) {
            if(mth_all.threads[i].state != MTH_STATE_ERROR) {
                pthread_cancel(mth_all.threads[i].thread);
            }
        }
        
        free((void *) mth_all.threads);
        mth_all.threads = NULL;
    }
    
    mth_all.nano = us_main * 1000;
    mth_all.count = threads;
    
    return ok;
}

bool mth_wait(void) {
    ulong nano = mth_all.nano;
    if(nano == 0) {
        return false;
    }
    
    struct timespec sleep = {.tv_sec = nano / 1000000000, .tv_nsec = nano % 1000000000};
    
    return nanosleep(&sleep, NULL) == 0;
}

static inline int mth_count(uint state) {
    if(! mth_lock()) {
        return -1;
    }
    
    int count = 0;
    for (int i = 0; i < mth_all.count; i++) {
        count += mth_all.threads[i].state == state ? 1 : 0;
    }
        
    return mth_unlock() ? count : -1;
}

uint mth_available(void) {
    return mth_count(MTH_STATE_WAIT);
}

bool mth_add_task(ulong ID, void *task) {
    ups_wait_online();
    
    if(! mth_lock()) {
        return false;
    }
    
    for (int i = 0; i < mth_all.count; i++) {
        if(mth_all.threads[i].state == MTH_STATE_WAIT) {
            mth_all.threads[i].ID = ID;
            mth_all.threads[i].task = task;
            mth_all.threads[i].state = MTH_STATE_READY;
                        
            return mth_unlock();
        }
    }
    
    mth_unlock();
       
    return false;
}

bool mth_batch(void **tasks, int len) {
    if(tasks == NULL || len <= 0 || mth_available() < mth_threads_count()) {
        return false;
    }
    
    for (int i = 0; i < len; i++) {
        if(tasks[i] == NULL) {
            return false;
        }
    }
    
    int started = 0, done = 0;
    
    bool ok = true;
    do {
        while(mth_results() > 0 && ok) {
            ok = ok && mth_get_result(NULL, NULL, NULL);
            done += ok ? 1 : 0;
        }
        
        while(started < len && mth_available() > 0) {
            ok = ok && mth_add_task(started, tasks[started]);
            started ++;
        }
        
        ok = ok && mth_errors() == 0;
        
        if(ok && done < len && mth_results() == 0) {
            mth_wait();
        }
    } while(done < len && ok);
        
    return ok && done == len; // redundant, not an issue
}

uint mth_results(void) {
    return mth_count(MTH_STATE_DONE);
}

bool mth_get_result(ulong *ID, void **result, void **task) {
    if(ID != NULL) {
        *ID = -1;
    }
    
    if(result != NULL) {
        *result = NULL;
    }
    
    if(task != NULL) {
        *task = NULL;
    }
    
    if(! mth_lock()) {
        return false;
    }
    
    for (int i = 0; i < mth_all.count; i++) {
        if(mth_all.threads[i].state == MTH_STATE_DONE) {
            if(ID != NULL) {
                *ID = mth_all.threads[i].ID;
            }
            
            if(result != NULL) {
                *result = (void *) mth_all.threads[i].result;
                mth_all.threads[i].result = NULL;
            }
            
            if(task != NULL) {
                *task = (void *) mth_all.threads[i].task;
                mth_all.threads[i].task = NULL;
            }
            
            mth_all.threads[i].state = MTH_STATE_WAIT;
            
            return mth_unlock();
        }
    }
    
    mth_unlock();
    
    return false;
}

bool mth_wait_results(bool all, bool clear, ulong usTimeOut) {
    ptime st, now;
    if(usTimeOut > 0) {
        clock_gettime(CLOCK_MONOTONIC, &st);
    }
    
    long us = 0;
    while (us <= usTimeOut) {
        int done = 0, wait = 0;
        
        if(! mth_lock()) {
            return false;
        }
        
        for (int i = 0; i < mth_all.count; i++) {
            byte state = mth_all.threads[i].state;
            
            if(state == MTH_STATE_ERROR) {
                mth_unlock();
                
                return false;
            }
            
            if(state == MTH_STATE_DONE) {
                if(mth_all.threads[i].result == NULL) {
                    mth_unlock();
                    
                    return false;
                }
                
                if(clear) {
                    mth_all.threads[i].result = NULL;
                    
                    mth_all.threads[i].state = MTH_STATE_WAIT;
                }
                
                if(! all) {
                    return mth_unlock();
                }
                
                done ++;
            } else if(state == MTH_STATE_WAIT) {
                wait ++;
            }
        }
        
        if(done + wait == mth_all.count) {
            return mth_unlock();
        }
        
        if(! mth_unlock()) {
            return false;
        }
        
        mth_wait();
        
        if(usTimeOut > 0) {
            clock_gettime(CLOCK_MONOTONIC, &now);
            
            us = now.tv_sec;
            us -= st.tv_sec;
            us *= MILLION;
            us += now.tv_nsec / 1000;
            us -= st.tv_nsec / 1000;
        }
    }
    
    return false;
}

uint mth_errors(void) {
    return mth_count(MTH_STATE_ERROR);
}

bool mth_get_error(ulong *ID, void **task) {
    if(! mth_lock()) {
        return false;
    }
    
    for (int i = 0; i < mth_all.count; i++) {
        if(mth_all.threads[i].state == MTH_STATE_ERROR) {
            if(ID != NULL) {
                *ID = mth_all.threads[i].ID;
            }
            
            if(task != NULL) {
                *task = (void *) mth_all.threads[i].task;
                mth_all.threads[i].task = NULL;
            }
            
            mth_all.threads[i].state = MTH_STATE_WAIT;
            
            return mth_unlock();
        }
    }
    
    mth_unlock();
    
    return false;
}

uint mth_running(void) {
    return mth_count(MTH_STATE_RUNNING);
}

bool mth_end_threads(bool ignore_errors, bool ignore_results) {
    if(mth_all.count == 0) { // there are no threads started
        return false;
    }
    
    uint count = mth_available() + (ignore_errors ? mth_errors() : 0) +
                (ignore_results ? mth_results() : 0);
    
    if(count < mth_all.count) { // there are uncollected results, unchecked errors or threads running
        return false;
    }
    
    mthread *th = (mthread *) mth_all.threads;
    for (int i = 0; i < mth_all.count; i++) {
        set_state(th + i, MTH_STATE_EXIT);
    }
    
    int ended;
    void *task;
    do {
        ended = 0;
        for (int i = 0; i < mth_all.count; i++) {
            if(get_state(th + i) == MTH_STATE_ENDED) {
                int status = pthread_join(th[i].thread, &task);
                set_state(th + i, MTH_STATE_JOINED);
                
                if(status != 0 || task != th + i) {
                    return false;
                }
                
                ended ++;
            } else if(get_state(th + i) == MTH_STATE_JOINED) {
                ended ++;
            }
        }
        
        if(ended < mth_all.count) {
            mth_wait();
        }
    } while(ended < mth_all.count);
    
    free((void *) mth_all.threads);
    mth_all.threads = NULL;
    mth_all.count = 0;
    mth_all.nano = 0;
    
    pthread_mutex_destroy(&mth_all.lock);
    
    return true;
}

bool mth_kill_threads(void) {
    if(mth_all.count == 0) {
        return false;
    }
    
    if(! mth_end_threads(true, true)) {
        for (int i = 0; i < mth_all.count; i++) {
            pthread_cancel(mth_all.threads[i].thread);
        }
        
        bool ok = true;
        for (int i = 0; i < mth_all.count; i++) {
            int status = pthread_join(mth_all.threads[i].thread, NULL);
            
            ok = ok && status == 0;
        }
        
        free((void *) mth_all.threads);
        mth_all.threads = NULL;
        mth_all.count = 0;
        mth_all.nano = 0;
        
        pthread_mutex_destroy(&mth_all.lock);
        
        return ok;
    }
    
    pthread_mutex_destroy(&mth_all.lock);
    
    return true;
}
