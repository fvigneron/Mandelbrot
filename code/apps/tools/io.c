//
//  io.c
//  Mandelbrot
//
//  Created by MIHALACHE Nicolae on 11/19/19.
//  Revised by VIGNERON François on 01/14/21.
//  Copyright © 2019 MIHALACHE Nicolae. All rights reserved.
//

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include <stdarg.h>
#include <mpfr.h>

#include "io.h"
#include "stopWatch.h"
#include "memFile.h"

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Pretty printing of numbers
// //////////////////////////////////////////////////////////////////////////////////////////
// TODO: what is the theoretical limitation for the size of c ?
// Formating the biggest long double takes 330 bytes, so 351 bytes with .20 decimals
// Pretty printing an fp80 could thus produce about 705 (not significant) bytes.
// If |c|>10^{76} then a premature truncations of the imaginary part may occur
// but attempting a pretty print on number that exceeds the precision of its
// storage structure is meaningless. Strings are allocated for a reasonable usage.
// TODO: why is the stu buffer shorter than the one for stc ? Update comment in io.h if change.

char *str(mpfr_t x) {
    char *s = malloc(200);
    mpfr_snprintf(s, 195, "%.50Rg", x);
    
    return s;
}

char *stl(fp80 c) {
    char *s = malloc(100);
    snprintf(s, 95, "%.20Lg %.20Lg", c->x, c->y);
    
    return s;
}

char *stc(mpc c) {
    char *s = malloc(200);
    mpfr_snprintf(s, 195, "%.40Rg %.40Rg", c->x, c->y);
    
    return s;
}

char *stv(mpv v, long i) {
    mpfr_t x;
    mpfr_init2(x, 150);
    
    mpv_get(x, v, i);
    
    char *res = str(x);
    mpfr_clear(x);
    
    return res;
}

char *stw(mpv v, long i) {
    mpc z;
    mpc_init(z, 150);
    
    mpv_getc(z, v, i);
    
    char *res = stc(z);
    mpc_clear(z);
    
    return res;
}

char *stbs(void *data, int len) {
    char *s = malloc(3 * len + 1);
    int sp = 0;
    byte *d = data;
    
    for (int i = 0; i < len; i++) {
        snprintf(s + sp, 3, "%02X", d[i]);
        sp += 2;
        
        if((i + 1) % 4 == 0) {
            if((i + 1) % 64 == 0) {
                snprintf(s + sp, 2, "\n");
                sp ++;
            } else if((i + 1) % 64 == 0) {
                snprintf(s + sp, 3, "  ");
                sp += 2;
            } else {
                snprintf(s + sp, 3, " ");
                sp ++;
            }
        }
    }
    
    return s;
}

char *sti(mpi x) {
    mpfr_t c, r;
    mpfr_init2(c, mpi_prec(x));
    mpfr_init2(r, 64);
    
    mpi_disk(c, r, x);
    
    char *s = malloc(200);
    mpfr_snprintf(s, 195, "%.40Rg %.10Rg", c, r);
    
    mpfr_clear(c);
    mpfr_clear(r);
    
    return s;
}

char *stu(u128 c) {
    mpfr_t x, y;
    mpfr_init2(x, 128);
    mpfr_init2(y, 128);
    
    u128_getr(x, y, c);
    char *s = malloc(200);
    mpfr_snprintf(s, 195, "%.40Rg %.40Rg", x, y);
    
    mpfr_clear(x);
    mpfr_clear(y);
    
    return s;
}

char *std(mpd d) {
    char *s = malloc(400);
    mpfr_snprintf(s, 395, "%.40Rg %.40Rg %.5Rg", d->x, d->y, d->r);
    
    return s;
}

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: fprintf wrappers
// //////////////////////////////////////////////////////////////////////////////////////////

static char error_file[IO_FILE_NAME_LEN] = "errors.txt";

bool io_error_file(char *file_name, bool rewrite) {
    if(file_name == NULL) {
        return false;
    }
    
    bool ok = snprintf(error_file, IO_FILE_NAME_LEN, "%s", file_name) < IO_FILE_NAME_LEN;
    if(ok && rewrite) {
        FILE *f = fopen(error_file, "w");
        if(f == NULL) {
            return false;
        }
        
        ok = ok && fclose(f) == 0;
    }
    
    return ok;
}

bool io_error(const char *format, ...) {
    va_list argp;
    
    if(error_file[0] == 0 || format == NULL) {
        return false;
    }
        
    FILE *f = fopen(error_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_error_echo(const char *format, ...) {
    va_list argp;
    
    va_start(argp, format);
    mpfr_vprintf(format, argp);
    va_end(argp);
    
    fflush(stdout);
    
    if(error_file[0] == 0 || format == NULL) {
        return false;
    }
    
    FILE *f = fopen(error_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

static char log_file[IO_FILE_NAME_LEN] = "log.txt";

bool io_log_file(char *file_name, bool rewrite) {
    if(file_name == NULL) {
        return false;
    }
    
    bool ok = snprintf(log_file, IO_FILE_NAME_LEN, "%s", file_name) < IO_FILE_NAME_LEN;
    if(ok && rewrite) {
        FILE *f = fopen(log_file, "w");
        if(f == NULL) {
            return false;
        }
        
        ok = ok && fclose(f) == 0;
    }
    
    return ok;
}

bool io_log(const char *format, ...) {
    va_list argp;
    
    if(log_file[0] == 0 || format == NULL) {
        return false;
    }
        
    FILE *f = fopen(log_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_log_echo(const char *format, ...) {
    va_list argp;
    
    va_start(argp, format);
    mpfr_vprintf(format, argp);
    va_end(argp);
    
    fflush(stdout);
    
    if(log_file[0] == 0 || format == NULL) {
        return false;
    }
    
    FILE *f = fopen(log_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

static char out_file[IO_FILE_NAME_LEN] = "out.txt";

bool io_out_file(char *file_name, bool rewrite) {
    if(file_name == NULL) {
        return false;
    }
    
    bool ok = snprintf(out_file, IO_FILE_NAME_LEN, "%s", file_name) < IO_FILE_NAME_LEN;
    if(ok && rewrite) {
        FILE *f = fopen(out_file, "w");
        if(f == NULL) {
            return false;
        }
        
        ok = ok && fclose(f) == 0;
    }
    
    return ok;
}

bool io_out(const char *format, ...) {
    va_list argp;
    
    if(out_file[0] == 0 || format == NULL) {
        return false;
    }
        
    FILE *f = fopen(out_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_out_echo(const char *format, ...) {
    va_list argp;
    
    va_start(argp, format);
    mpfr_vprintf(format, argp);
    va_end(argp);
    
    fflush(stdout);
    
    if(out_file[0] == 0 || format == NULL) {
        return false;
    }
    
    FILE *f = fopen(out_file, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_write(const char *fileName, const char *format, ...) {
    va_list argp;
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_append(const char *fileName, const char *format, ...) {
    va_list argp;
    
    FILE *f = fopen(fileName, "a");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool io_fprint(const char *fileName, bool append, bool echo, const char *format, ...) {
    va_list argp;
    
    if(echo) {
        va_start(argp, format);
        mpfr_vprintf(format, argp);
        va_end(argp);
        
        fflush(stdout);
    }
    
    if(fileName == NULL) {        
        return false;
    }
    
    FILE *f = fopen(fileName, append ? "a" : "w");
    if(f == NULL) {
        return false;
    }
    
    va_start(argp, format);
    int count = mpfr_vfprintf(f, format, argp);
    va_end(argp);
    
    int close = fclose(f);
    
    return count > 0 && close == 0;
}

bool dir(char dirName[]) {
    if(access(dirName, F_OK) != 0) {
        if(mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO) != 0) {
            printf("Could not create dir %s\n", dirName);
            
            return false;
        }
    }
    
    return true;
}

bool file_can_read(char *fn) {
    return access(fn, R_OK) != 0;
}

bool file_can_write(char *fn) {
    return access(fn, W_OK) != 0;
}

bool file_exists(char *fn) {
    return access(fn, F_OK) != 0;
}


long io_file_size(FILE *f) {
    if(f == NULL) {
        return 0;
    }
    
    long pos = ftell(f);
    if(pos < 0) {
        return 0;
    }

    // check the size of the file
    bool ok = fseek(f, 0, SEEK_END) == 0;
    long fs = ftell(f);
    ok = ok && fseek(f, pos, SEEK_SET) == 0;

    return ok ? fs : -1;
}

bool io_file_seek(FILE *f, long pos, bool extend) {
    if(f == NULL || pos < 0) {
        return false;
    }
    
    long size = io_file_size(f);
    if(size >= pos) {
        bool ok = fseek(f, pos, SEEK_SET) == 0;
        ok = ok && ftell(f) == pos;
        
        return ok;
    }
    
    if(extend && ftruncate(fileno(f), pos) != 0) {
        return false;
    }
    
    return io_file_seek(f, pos, false);
}

long io_size_of(char *fileName) {
    FILE *f = fopen(fileName, "r");
    
    if(f == NULL) {
        return -1;
    }
    
    long fs = io_file_size(f);
    fclose(f);
    
    return fs;
}

long io_file_remainder(FILE *f) {
    if(f == NULL) {
        return 0;
    }
    
    long pos = ftell(f);
    if(pos < 0) {
        return 0;
    }

    // check the size of the file
    bool ok = fseek(f, 0, SEEK_END) == 0;
    long fs = ftell(f);
    ok = ok && fseek(f, pos, SEEK_SET) == 0;

    return ok && fs >= pos ? fs - pos : -1;
}

void static fmt(double v, char *str, char *unit) {
    int d = v < 9.995 ? 2 : v < 99.95 ? 1 : 0;
    snprintf(str, 19, "%%.%dlf %s", d, unit);
}

void mem_size(ulong mem, char *str) {
    char num[20];
    double x;
    
    if(mem < 1 << 10) {
        snprintf(str, 29, "%ld bytes", mem);
    } else if(mem < 1 << 20) {
        x = ldexp(mem, -10);
        fmt(x, num, "KB");
        snprintf(str, 29, num, x);
    } else if(mem < 1 << 30) {
        x = ldexp(mem, -20);
        fmt(x, num, "MB");
        snprintf(str, 29, num, x);
    } else if(mem < 1L << 40) {
        x = ldexp(mem, -30);
        fmt(x, num, "GB");
        snprintf(str, 29, num, x);
    } else {
        x = ldexp(mem, -40);
        fmt(x, num, "TB");
        snprintf(str, 29, num, x);
    }
}

void file_size(ulong len, char *str) {
    if(! FILE_SIZE_DECIMAL) {
        mem_size(len, str);
        
        return;
    }
    
    char num[20];
    double x;
    
    if(len < 1 << 10) {
        snprintf(str, 29, "%ld bytes", len);
    } else if(len < 1 << 20) {
        x = len * 1E-3;
        fmt(x, num, "KB");
        snprintf(str, 29, num, x);
    } else if(len < 1 << 30) {
        x = len * 1E-6;
        fmt(x, num, "MB");
        snprintf(str, 29, num, x);
    } else if(len < 1L << 40) {
        x = len * 1E-9;
        fmt(x, num, "GB");
        snprintf(str, 29, num, x);
    } else {
        x = len * 1E-12;
        fmt(x, num, "TB");
        snprintf(str, 29, num, x);
    }
}

bool disk_stats(char *path) {
    struct statvfs fiData;
    
    if((statvfs(path, &fiData)) < 0 ) {
        printf("Failed to stat %s:\n", path);
        
        return false;
    } else {
        printf("Disk %s: \n", path);
        printf("\tblock size: %lu\n", fiData.f_frsize);
        printf("\ttotal no blocks: %lu\n", (ulong) fiData.f_blocks);
        printf("\tfree blocks: %lu\n", (ulong) fiData.f_bfree);
        
        return true;
    }
}

long disk_free(char *path) {
    struct statvfs fiData;
    
    if((statvfs(path, &fiData)) < 0 ) {
        return -1;
    } else {
        long free = fiData.f_frsize;
        free *= fiData.f_bfree;
        
        return free;
    }
}

size_t fread_block(void *ptr, size_t size, size_t nmemb, size_t block, FILE *stream) {
    size_t read = 0;
    size_t step = block / size;
    if(step <= 0 || ptr == NULL || stream == NULL) {
        return 0;
    }
    
    while(read + step <= nmemb) {
        size_t lr = fread(ptr + read * size, size, step, stream);
        read += lr;
        
        if(lr < step) {
            return read;
        }
    }
    
    if(read < nmemb) {
        read += fread(ptr + read * size, size, nmemb - read, stream);
    }
    
    return read;
}

void io_print_start(void) {
    char now[100];
    date(now, 100);
    
    printf("Started on %s\n", now);
}

bool io_update_md5(MD5_CTX *md5, FILE *f, long size) {
    if(md5 == NULL || f == NULL || size <= 0) {
        return false;
    }
    
    long mem = size < FILE_BLOCK ? size : FILE_BLOCK;
    byte *buf = malloc(mem);
    if(buf == NULL) {
        return false;
    }
    
    bool ok = true;
    long done = 0;
    
    while(ok && done < size) {
        long read = size - done > mem ? mem : size - done;
        ok = ok && fread(buf, 1, read, f) == read;
        ok = ok && MD5_Update(md5, buf, (uint) read);
        
        done += read;
    }
    
    free(buf);
    
    return ok;
}

static int csv_line_count(char *line) {
    if(line == NULL || line[0] == '#' || line[0] == ';' || line[0] == '/') {
        return 0;
    }
    
    int cols = 0, i = 0;
    bool inq = false;
    while(line[i] != 0 && i < CSV_LINE_MAX_LEN) {
        char c = line[i++];
        
        if(inq) {
            if(c != '"') {
                continue;
            }
            
            inq = false;
        }
        
        if(c == ',') {
            if(cols == 0) {
                cols = 2;
            } else {
                cols ++;
            }
        }
        
        if(cols == 0 && c != ' ') {
            cols = 1;
        }
        
        if(c == '"') {
            inq = true;
        }
    }
    
    return cols;
}

static void trim(strings str, int index) {
    int pos = 0;
    char *line = str->fields[index];
    
    // trim the beginning
    while(line[pos] == ' ') {
        pos ++;
    }
    str->fields[index] = str->fields[index] + pos;
    int st = pos;
    
    // trim the end
    while(line[pos] != 0) {
        pos ++;
    }
    pos --;
    while(pos >= 0 && (line[pos] == ' ' || line[pos] == '\n')) {
        pos --;
    }
    line[pos + 1] = 0;
    
    // save the count
    str->lens[index] = pos - st;
}

static bool csv_line_split(strings l) {
    if(l == NULL || l->count <= 0) {
        return false;
    }
    
    l->fields[0] = l->line;
    
    int cols = 1, i = 0;
    bool inq = false;
    while(cols <= l->count && l->line[i] != 0 && i < CSV_LINE_MAX_LEN) {
        char c = l->line[i++];
        
        if(inq) {
            if(c != '"') {
                continue;
            }
            
            l->line[i - 1] = ' '; // replace quotes by spaces
            inq = false;
        }
        
        if(c == ',') {
            l->fields[cols] = l->line + i;
            l->line[i - 1] = 0;   // end the previous field
            cols ++;
            
            // trim spaces from the previous column
            trim(l, cols - 2);
        }
        
        if(c == '"') {
            l->line[i - 1] = ' '; // replace quotes by spaces
            inq = true;
        }
    }
    
    // trim spaces from the last column
    trim(l, cols - 1);
    
    return cols == l->count && l->line[i] == 0;
}

strings io_read_csv_line(FILE *f) {
    if(f == NULL || feof(f)) {
        return NULL;
    }
    
    strings line = malloc(sizeof(strings_struct));
    line->line = malloc(sizeof(char) * CSV_LINE_MAX_LEN);
    line->fields = NULL;
        
    do {
        char *l = fgets(line->line, CSV_LINE_MAX_LEN, f);
        
        if(l != line->line) {
            strings_free(line);
            
            return NULL;
        }
        
        line->count = csv_line_count(line->line);
    } while(line->count == 0);
    
    line->fields = malloc(line->count * sizeof(char *));
    line->lens = malloc(line->count * 4);
    
    if(! csv_line_split(line)) {
        strings_free(line);
        
        return NULL;
    }
    
    return line;
}

void strings_free(strings line) {
    if(line == NULL) {
        return;
    }
    
    if(line->fields != NULL) {
        free(line->fields);
        line->fields = NULL;
    }
    
    if(line->lens != NULL) {
        free(line->lens);
        line->lens = NULL;
    }
    
    if(line->line != NULL) {
        free(line->line);
        line->line = NULL;
    }
    
    line->count = 0;
    
    free(line);
}

bool str_starts_with(char *str, char *start) {
    if(str == NULL || start == NULL) {
        return false;
    }
    
    ulong ll = strlen(str);
    ulong sl = strlen(start);
    
    return sl <= ll && strncmp(str, start, sl) == 0;
}

bool strn_starts_with(char *str, int ll, char *start, int sl) {
    if(sl == 0) {
        return true;
    }
    
    if(str == NULL || start == NULL) {
        return false;
    }
    
    return sl <= ll && strncmp(str, start, sl) == 0;
}

bool str_ends_with(char *str, char *end) {
    if(str == NULL || end == NULL) {
        return false;
    }
    
    ulong ll = strlen(str);
    ulong sl = strlen(end);
    
    return sl <= ll && strncmp(str + (ll - sl), end, sl) == 0;
}

bool strn_ends_with(char *str, int ll, char *end, int sl) {
    if(sl == 0) {
        return true;
    }
    
    if(str == NULL || end == NULL) {
        return false;
    }
    
    return sl <= ll && strncmp(str + (ll - sl), end, sl) == 0;
}


strings io_list_dir(char *dir, bool recursive) {
    return io_filter_dir(dir, NULL, NULL, recursive);
}

bool io_add_string_to_mfile(mfile f, char *str, int len) {
    bool ok = mfile_puti(f, len);
    ok = ok && mfile_putbs(f, (byte *) str, len + 1);
    
    return ok;
}

strings io_mfile_to_strings(mfile f, int count) {
    if(f == NULL) {
        return NULL;
    }
    
    if(count == 0) {
        mfile_free(f);
        
        return NULL;
    }
    
    strings list = malloc(sizeof(strings_struct));
    list->line = malloc(f->pos - (count << 2));
    list->fields = malloc(sizeof(char *) * count);
    list->lens = malloc(4 * count);
    
    list->count = count;
    f->pos = 0;
    
    ulong pos = 0;
    bool ok = true;
    for (int i = 0; i < count && ok; i++) {
        uint len = mfile_geti(f);
        list->lens[i] = len;
        list->fields[i] = list->line + pos;
        ok = ok && mfile_getbs(f, (byte *) list->line + pos, len + 1);
        pos += len + 1;
    }
    
    ok = ok && pos == f->len - (count << 2);
    mfile_free(f);
    
    if(! ok) {
        strings_free(list);
        list = NULL;
    }
    
    return list;
}

strings io_filter(strings list, char *start, char *end) {
    if(list == NULL) {
        return NULL;
    }
    
    mfile f = mfile_new(1000);
    int count = 0;
    bool ok = true;
    
    int sl = start == NULL ? 0 : (int) strlen(start);
    int el = end == NULL ? 0 : (int) strlen(end);
    
    for (int i = 0; i < list->count && ok; i++) {
        if(! strn_starts_with(list->fields[i], list->lens[i], start, sl) ||
           ! strn_ends_with(list->fields[i], list->lens[i], end, el)) {
            continue;
        }
        
        ok = ok && io_add_string_to_mfile(f, list->fields[i], list->lens[i]);
        count ++;
    }
    
    if(! ok) {
        mfile_free(f);
        
        return NULL;
    }
    
    return io_mfile_to_strings(f, count);
}

strings io_filter_dir(char *dir, char *start, char *end, bool recursive) {
    if(dir == NULL) {
        return NULL;
    }
    
    mfile f = mfile_new(1000);
    DIR *d = opendir(dir);
    if(d == NULL) {
        return NULL;
    }
    
    int sl = start == NULL ? 0 : (int) strlen(start);
    int el = end == NULL ? 0 : (int) strlen(end);
    
    int count = 0;
    bool ok = true;
    struct dirent *dirt;
    while ((dirt = readdir(d)) != NULL && ok) {
        char *s = dirt->d_name;
        
        if(recursive && dirt->d_type == DT_DIR && strcmp(s, ".") != 0 && strcmp(s, "..") != 0) {
            ulong n = strlen(dir) + strlen(s) + 2;
            char ndir[n];
            snprintf(ndir, n, "%s%s" FOLDER_SEPARATOR, dir, s);
            
            strings chld = io_filter_dir(ndir, start, end, true);
            if(chld != NULL) {
                for (int i = 0; i < chld->count; i++) {
                    ok = ok && io_add_string_to_mfile(f, chld->fields[i], chld->lens[i]);
                    count ++;
                }
                
                strings_free(chld);
            }
        }
        
        int ll = (int) strlen(s);
        
        if(! strn_starts_with(s, ll, start, sl) || ! strn_ends_with(s, ll, end, el)) {
            continue;
        }
        
        ulong n = strlen(dir) + strlen(s) + 2;
        char fn[n];
        snprintf(fn, n, "%s%s", dir, s);
        
        ok = ok && io_add_string_to_mfile(f, fn, (uint) strlen(fn));
        count ++;
    }
    closedir(d);
    
    if(! ok) {
        mfile_free(f);
        
        return NULL;
    }
    
    return io_mfile_to_strings(f, count);
}

void io_print_strings(strings list) {
    if(list == NULL) {
        return;
    }
    
    for (int i = 0; i < list->count; i++) {
        printf("%s\n", list->fields[i]);
    }
}

bool io_cat(char *dst, char **src, int count, bool delete_src) {
    if(dst == NULL || src == NULL || count <= 0) {
        return false;
    }
    
    bool ok = true;
    for (int i = 0; i < count && ok; i++) {
        mfile f = mfile_load(src[i]);
        ok = ok && f != NULL;
        ok = ok && mfile_save(dst, f, i > 0);
        
        if(delete_src) {
            ok = ok && remove(src[i]) == 0;
        }
        
        mfile_free(f);
    }
    
    return ok;
}
