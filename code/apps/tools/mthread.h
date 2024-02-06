//
//  mthread.h
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

/**
 \file mthread.h
 \brief Data structure and methods for multi-threading with a task consumer model.
 
 This is based on "pthread.h" and offers a simplification of thread-management tasks. All methods
 in this module should be called from the main thread.
 
 There is ony one instance of this class / module, so for different types of tasks that need to
 distribute tasks to several computing threads, @ref mth_init() and @ref mth_end_threads()
 (or @ref mth_kill_threads()) have to be used sequentially.
 
 While waiting for the computing tasks to complete, the main thread should use @ref mth_wait() to reduce
 CPU core utilisation.
*/
#ifndef mthread_h
#define mthread_h

#include "ntypes.h"

// MARK: Public types

/// A method that computes a @c task and returns a result.
typedef void *(* mth_run)(void *task, long task_ID, int thread_index);

// MARK: Private definitions, do not use directly

/// The thread will sleep and check periodically until the state becomes @c MTH_STATE_READY or @c MTH_STATE_EXIT
#define MTH_STATE_WAIT     0

/// The task to be done is prepared, ready to be executed
#define MTH_STATE_READY    1

/// The thread is executing the task
#define MTH_STATE_RUNNING  2

/// The task has been completed successfully
#define MTH_STATE_DONE     3

/// Some error occurred and the taks did not complete successfully
#define MTH_STATE_ERROR    4

/// All jobs done, the thread should end
#define MTH_STATE_EXIT     5

/// All jobs done, the thread should end
#define MTH_STATE_ENDED    6

/// The thread has physically ended, was joined. Do NOT attempt to join again !
#define MTH_STATE_JOINED   7

/// The maximum number of threads
#define MTH_MAX_THREADS   72

/// A key to ensure memory synch at the start of a thread.
#define MTH_KEY  ((ulong) ((1L << 62) * PI))

typedef struct {
    volatile byte state;          ///< the state of the thread, one of @c MTH_STATE_XXX
    volatile int index;           ///< the index of the thread
    volatile long ID;             ///< optional ID provided by the user
    volatile void *task;          ///< data needed for perfroming the task
    volatile void *result;        ///< the result of the task
    volatile ulong nano;          ///< nanoseconds to sleep between tasks
    volatile mth_run run;         ///< the specific run function
    volatile pthread_t thread;    ///< the thread
    volatile ulong key;           ///< a key to insure the proper start of a thread after memory synch
} mthread;

typedef struct {
    volatile uint count;                ///< the number of available threads
    volatile ulong nano;                ///< nanoseconds to sleep when asked to wait
    volatile mthread *threads;          ///< the list of tasks
    pthread_mutex_t lock;               ///< the synchronization lock
} mth_threads;

// MARK: Public methods

/// @brief Initializes and starts @c threads threads that will execute the method @c run repeatedly.
///
/// While @ref mth_available() > 0, @ref mth_add_task() could be used to submit a new computing task.
/// There are new results while @ref mth_results() > 0 and in this case they can be retrieved by @ref mth_get_result().
/// End all tasks with @ref mth_end_threads() or with @ref mth_kill_threads().
///
/// @param threads the number of threads to start
/// @param run the method that computes one task and exits afterwards
/// @param micros_thread the time in microseconds for threads to sleep, while waiting for new computing tasks
/// @param micros_main the time in microseconds for the main thread to sleep, when invoking @ref mth_wait()
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_init(uint threads, mth_run run, ulong micros_thread, ulong micros_main);

/// @brief Returns the number of threads that have been started.
///
/// @return the number of threads that have been started
uint mth_threads_count(void);

/// @brief Sleeps for the duration specified at init by @ref mth_init().
///
/// @return @ref true if the thread was suspended for the entire duration, @ref false otherwise
bool mth_wait(void);

/// @brief Returns the number of threads that are waiting for a computing task.
///
/// @return the number of threads that are waiting for a computing task
uint mth_available(void);

/// @brief Distrubutes the computing @c task to the first available thread.
///
/// @warning Fails if all threads are running, have uncollected results or unchecked errors.
/// Check if @ref mth_available() > 0 before calling this method.
///
/// @param ID the id of the computing task
/// @param task the data necessary to compute the task
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_add_task(ulong ID, void *task);

/// @brief Distrubutes the list of computing @c tasks to all threads.
///
/// @warning This method is blocking until an error has occurred or all tasks in the list have been processed.
/// @warning This method fails if there are other threads running, uncollected results or unchecked errors.
///
/// @param tasks the list of tasks
/// @param len the length of the list
///
/// @return @ref true if all tasks have completed sucessfully, @ref false otherwise
bool mth_batch(void **tasks, int len);

/// @brief Returns the number of threads that have completed the computing task.
///
/// @return the number of threads that have completed the computing task
uint mth_results(void);

/// @brief Retrieves the first available result and stores it to @c result.
///
/// Fails if there are no available results. Check if @ref mth_results() > 0 before calling this method.
/// Does not fail if any of the provided pointers is  @c NULL, ignore it on a case by case basis.
///
/// @param ID the pointer where to retrieve the id of the computing task
/// @param result the pointer where to retrieve the result
/// @param task the pointer where to retrieve the the data provided to compute the task
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_get_result(ulong *ID, void **result, void **task);

/// @brief Waits in the main thread for results. Fails if there are errors or if the @c usTimeOut is exceeded.
///
/// If @c usTimeOut==0, there is no time limit. A completed task that has returned a @c NULL result is
/// considered to be an error.
///
/// @param all @ref true to wait for all results, @ref false to return as soon as there is one available
/// @param clear @ref true to clear available results, @ref false to preserve them
/// @param usTimeOut the time limit in micros-seconds
///
/// @return @ref true if results are available, @ref false if some error occurred or the time limit is exceeded.
bool mth_wait_results(bool all, bool clear, ulong usTimeOut);

/// @brief Returns the number of threads that could not complete the computing task.
///
/// @return the number of threads that could not complete the computing task
uint mth_errors(void);

/// @brief Retrieves the first available error and stores its ID to @c ID.
///
/// Fails if there are no errors. Check if @ref mth_errors() > 0 before calling this method.
/// Does not fail if any of the provided pointers is  @c NULL, ignore it on a case by case basis.
///
/// @param ID the pointer where to retrieve the id of the computing task
/// @param task the pointer where to retrieve the the data provided to compute the task
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_get_error(ulong *ID, void **task);

/// @brief Returns the number of threads that are processing a computing task.
///
/// @return the number of threads that are processing a computing task
uint mth_running(void);

/// @brief If there are no runing computing tasks, signals to all threads that they should exit.
///
/// Fails if there are running computing tasks or if there are results or errors that are not ignored.
///
/// @param ignore_errors @true to ignore existing errors and proceed with the cleaup, @ref false orherwise
/// @param ignore_results @true to ignore available results and proceed with the cleaup, @ref false orherwise
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_end_threads(bool ignore_errors, bool ignore_results);

/// @brief Sends a hard system signal to all threads to stop executing.
///
/// Fails only if there are no initialized threads, but does not check if the threads have effectively stopped running.
///
/// @return @ref true if sucessfull, @ref false otherwise
bool mth_kill_threads(void);

#endif /* mthread_h */
