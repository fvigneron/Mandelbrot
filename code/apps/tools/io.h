//
//  io.h
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

 /**
  \file io.h
  \brief A collection of basic @c IO functions, mostly used for debugging.
 */

#ifndef io_h
#define io_h

#include <mpfr.h>
#include <sys/statvfs.h>

#include "ntypes.h"
#include "mpc.h"
#include "mpi.h"
#include "mpd.h"
#include "mpv.h"
#include "u128c.h"


/// The default upper bound of the size of a block read from a file in one system call.
#define FILE_BLOCK (1L << 26)

/// Base in which kB, MB, GB and TB are displayed for files.
#define FILE_SIZE_DECIMAL 0

/// Base in which kB, MB, GB and TB are displayed for files.
#define CSV_LINE_MAX_LEN 100000

#define IO_FILE_NAME_LEN 250

#ifdef _WIN32
#define FOLDER_SEPARATOR "\\"
#else
#define FOLDER_SEPARATOR "/"
#endif

/// @struct csv_line_struct
/// @brief A line in a CSV file.
typedef struct {
    int count;         ///< the number of columns
    char *line;        ///< the line, with extra '0' characters to delimit fields
    char **fields;     ///< the fileds, as pointers to substrings of @c line
    int *lens;         ///< the legths of the fields
} strings_struct;

/// Convenience pointer to csv_line_struct
typedef strings_struct *strings;

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: Pretty printing of numbers
// //////////////////////////////////////////////////////////////////////////////////////////

/// @brief Pretty print an @c mpfr_t real number into a string.
///
/// Writes into a string the real number @c x, with 50 decimals after the decimal point.
///
/// Returns a pointer to that string. The buffer string is 200 bytes long.
/// @param x @c mpfr_t real number
char *str(mpfr_t x);

/// @brief Pretty print an @c fp80 complex number into a string.
///
/// Writes into a string the real and imaginary part of @c c, with 20 decimals after the decimal point and with a space as separator.
///
/// Returns a pointer to that string. The buffer string is 100 bytes long.
/// @param c @c fp80 complex number
char *stl(fp80 c);

/// @brief Pretty print an @c mpc complex number into a string.
///
/// Writes into a string the real and imaginary part of @c c, with 40 decimals after the decimal point and with a space as separator.
///
/// Returns a pointer to that string. The buffer string is 200 bytes long.
/// @param c @c mpc complex number
char *stc(mpc c);

/// @brief Pretty print the interval @c c into a string.
///
/// Writes into a string the center and radius of the interval @c c, with 40 and respectively 10 decimals after the decimal point
/// and with a space as separator.
///
/// Returns a pointer to that string. The buffer string is 200 bytes long.
/// 
/// @param c @c mpc complex number
char *sti(mpi c);

/// @brief Pretty print an @c u128 complex number into a string.
///
/// Writes into a string the real and imaginary part of @c c, with 40 decimals after the decimal point and with a space as separator.
///
/// Returns a pointer to that string. The buffer string is 200 bytes long.
/// @param c @c u128 complex number
char *stu(u128 c);

/// @brief Pretty print an @c mpd complex disk into a string.
///
/// Writes into a string the real and imaginary part of the center of @c d, with 40 decimals after the decimal point,
/// followed by the radius of @c d with 5 decimals after the decimal point and with a space as separator.
///
/// Returns a pointer to that string. The buffer string is 400 bytes long.
/// @param d @c mpd complex disk
char *std(mpd d);

/// @brief Pretty print an element of a @c mpv vector into a string.
///
/// Writes into a string the real number @c v[i], with 50 decimals after the decimal point.
///
/// @param v @c mpv vector
/// @param i the position in the vector
///
/// @return a pointer to that string. The buffer string is 400 bytes long.
char *stv(mpv v, long i);

/// @brief Pretty print an element of a @c mpv vector into a string.
///
/// Writes into a string the complex number @c v[i], with 40 decimals after the decimal point for each coordinate.
///
/// @param v @c mpv vector
/// @param i the position in the vector
///
/// @return a pointer to that string. The buffer string is 400 bytes long.
char *stw(mpv v, long i);

/// @brief Pretty hex print for a vector of bytes into a string.
///
/// @param data the data
/// @param len the number of bytes to convert
///
/// @return a pointer to that string, with length depending on @c len
char *stbs(void *data, int len);

// //////////////////////////////////////////////////////////////////////////////////////////
// MARK: fprintf wrappers
// //////////////////////////////////////////////////////////////////////////////////////////

bool io_error_file(char *file_name, bool rewrite);
bool io_error(const char *format, ...);
bool io_error_echo(const char *format, ...);

bool io_log_file(char *file_name, bool rewrite);
bool io_log(const char *format, ...);
bool io_log_echo(const char *format, ...);

bool io_out_file(char *file_name, bool rewrite);
bool io_out(const char *format, ...);
bool io_out_echo(const char *format, ...);

/// @brief Wrapper for @c mpfr_fprintf that writes into a fresh named file.
///
/// Writes the C string pointed to by @c format into the file pointed to by @c fileName.
/// Additional arguments following @c format are formatted and inserted in the usual way of @c mpfr_fprintf.
///
/// Automatically handles the opening and closing of the file @c fileName.
/// The file will either be created or overwritten (@c write mode).
///
/// Returns 1 if all operations succeded, 0 otherwise.
/// @param fileName String containing the name of the file to be used.
/// @param format String (with optional format specifiers) to be written in @c fileName.
///
/// @return @ref true if successfull, @ref false otherwise
bool io_write(const char *fileName, const char *format, ...);

/// @brief Wrapper for @c fprintf that writes at the end of a named file.
///
/// Writes the C string pointed to by @c format into the file pointed to by @c fileName.
/// Additional arguments following @c format are formatted and inserted in the usual way of @c fprintf.
///
/// Automatically handles the opening and closing of the file @c fileName.
/// The output will be written at the end of the file, which is created if necessary (@c append mode).
///
/// @param fileName String containing the name of the file to be used.
/// @param format String (with optional format specifiers) to be written in @c fileName.
///
/// @return @ref true if successfull, @ref false otherwise
bool io_append(const char *fileName, const char *format, ...);

/// @brief Wrapper for @c mpfr_fprintf that writes to a file, with an optional copy to @c stdout.
///
/// Writes the C string pointed to by @c format into the file pointed to by @c fileName.
/// Additional arguments following @c format are formatted and inserted in the usual way of @c mpfr_fprintf.
///
/// Automatically handles the opening and closing of the file @c fileName.
/// The output will be written at the end of the file if @c append, which is created if necessary (@c append mode),
/// the content is erased otherwise.
///
/// If @c echo, a copy of the message will be written to the standard output @c stdout.
///
/// @param fileName string containing the name of the file to be used.
/// @param append @ref true to append to the text file, @ref false to erase its content
/// @param echo toggles writing a copy to @c stdout on (if @c print) or off (if @c !print)
/// @param format string (with optional format specifiers) to be written in @c fileName
///
/// @return @ref true if successfull, @ref false otherwise
bool io_fprint(const char *fileName, bool append, bool echo, const char *format, ...);

/// @brief Created the folder with name @c dirName if it does not exist.
///
/// @param dirName the path of the folder
///
/// @return @ref true if the folder is ready to use , @ref false otherwise
bool dir(char dirName[]);

/// @brief Returns the size of the file, in bytes.
///
/// If the return value is positive or @c 0, it leaves the file a the same position.
///
/// @param f the file
///
/// @return the file size in bytes, @c -1 if the file does not exist or some error occurred
long io_file_size(FILE *f);

/// @brief Positions the file pointer in the file @c f to absolute position @c pos, while increasing the size of the file, if needed.
///
/// @param f the file
/// @param pos the desired absolute position in @c f
/// @param extend @ref true to extend the file if it is shorter than @c pos
///
/// @return @ref true if the file pointer is at the position @c pos, @ref false otherwise
bool io_file_seek(FILE *f, long pos, bool extend);

bool file_can_read(char *fn);
bool file_can_write(char *fn);
bool file_exists(char *fn);

/// @brief Returns the size of the file, in bytes.
///
/// @param fileName the path of the file
///
/// @return the file size in bytes, @c -1 if the file does not exist or some error occurred
long io_size_of(char *fileName);

/// @brief Returns the size of the file after the current position, in bytes.
///
/// If the return value is positive or @c 0, it leaves the file a the same position.
///
/// @param f the file
///
/// @return the relative file size from the current position,  in bytes, @c -1 if the file does not exist or some error occurred
long io_file_remainder(FILE *f);

/// @brief Pretty prints the memory size into the string @c str.
///
/// For memory sizes one KB is 1024 bytes, 1 MB = 2^20 bytes etc.
///
/// @warning The string should be at least @c 30 bytes long.
///
/// @param mem the memory size in bytes
/// @param str human readable string that will be stamped
void mem_size(ulong mem, char *str);

/// @brief Pretty prints the memory size into the string @c str.
///
/// For file sizes one KB is 1000 bytes, 1 MB = 1000000 bytes etc.
///
/// @warning The string should be at least @c 30 bytes long.
///
/// @param mem the memory size in bytes
/// @param str human readable string that will be stamped
void file_size(ulong mem, char *str);

/// @brief Prints stats of the volume containing the given @c path.
///
/// @param path a path on the desired volume
///
/// @return @ref true if successfull, @ref false otherwise
bool disk_stats(char *path);

/// @brief Returns the available space in bytes of the volume containing the given @c path.
///
/// @param path a path on the desired volume
///
/// @return the available space in bytes or @c -1 if the stat cannot be obtained
long disk_free(char *path);

/// @brief Replacement for the standard function @c fread to limit the size of file blocks read at one system call.
///
/// @param ptr the destination of the data read from the file
/// @param size the size of a record
/// @param nmemb the number of records to read
/// @param block the maximum size of a block read from the file with one system call
/// @param stream the input file
size_t fread_block(void *ptr, size_t size, size_t nmemb, size_t block, FILE *stream);

/// @brief Prints "Started on #date" with the current date and time.
void io_print_start(void);

bool io_update_md5(MD5_CTX *md5, FILE *f, long size);

/// @brief Reads the next line from @c f and splits it in the standard way.
///
/// If the line starts with '#', ';' or with '/', it is considered a comment and the next line is read.
/// The same is true if it contains only spaces.
///
/// @param f the input text file.
///
/// @return the @c csv_line or @c NULL if EOF is reached or some error occurred.
strings io_read_csv_line(FILE *f);

/// @brief Frees the memory used by a @c csv_line.
///
/// If @c line is @c NULL, it does not fail.
///
/// @param line the line.
void strings_free(strings line);

bool str_starts_with(char *str, char *start);
bool strn_starts_with(char *str, int len_str, char *start, int len_start);

bool str_ends_with(char *str, char *end);
bool strn_ends_with(char *str, int len_str, char *end, int len_end);

strings io_list_dir(char *dir, bool recursive);
strings io_filter_dir(char *dir, char *start, char *end, bool recursive);

bool io_add_string_to_mfile(mfile f, char *str, int len);
strings io_mfile_to_strings(mfile f, int count);
void io_print_strings(strings list);

bool io_cat(char *dst, char **src, int count, bool delete_src);

#endif /* io_h */
