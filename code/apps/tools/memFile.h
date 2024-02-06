//
//  memFile.h
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
 \file memFile.h
 \brief A memory buffer for a file that mimics the system functions fread(), fwrite() and fseek().
 
  It offers other convenience access functions and direct access to the data, to use in the inner loops of intensive IO operations.
*/

#ifndef memFile_h
#define memFile_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ntypes.h"

#ifdef __APPLE__
#  include <CommonCrypto/CommonDigest.h>
#  define MD5_DIGEST_LENGTH           CC_MD5_DIGEST_LENGTH
#  define MD5_CTX                     CC_MD5_CTX
#  define MD5_Init                    CC_MD5_Init
#  define MD5_Update                  CC_MD5_Update
#  define MD5_Final                   CC_MD5_Final
#else
#  include <openssl/md5.h>
#endif


// the minimal size of the header of a binary file
#define HEADER_MIN_LEN 20

// the number of bytes used to store the MD5 checksum in the header of a binary file
#define HEADER_MD5_LEN 16

/// Memory buffer increase ratio: doubles the capacity every two steps
#define MFILE_FACTOR  1.42
/// Memory buffer additive increase step, in bytes
#define MFILE_ADD     100

/// @struct memFile_struct
/// @brief The data structure for a memory file @c mfile.
typedef struct {
    ulong cap;  ///< the total allocated capacity, in bytes
    ulong len;  ///< the length of the file, at most @c cap
    ulong pos;  ///< the position in the file
    byte *data; ///< the data buffer
} memFile_struct;

/// @union fileHeaderId
/// \brief Data structure used for file type identification.
typedef union {
    char type[9];
    ulong fileTypeId;
} fileHeaderId;

/// Convenience type for a pointer to a memory file memFile_struct.
typedef memFile_struct *mfile;

/// The amount of bytes available to read from the file.
#define mfile_space(m) (m->len - m->pos)

/// Erases the content of the file.
#define mfile_erase(m) (m->len = 0)

/// Tests the end of file status.
#define meof(m) (m->len <= m->pos)

/// @brief Creates a new file with the given capacity @c cap.
///
/// @param cap the initial capacity of the memory buffer
///
/// @return the new file
mfile mfile_new(ulong cap);

/// @brief Frees all the memory used by the file @c m, assuming the struct has been allocated with @c malloc(), for example
/// with @c mfile_new().
///
/// It does @b not fail if @c ps==NULL.
///
/// @param m the file
///
/// @return @ref true if the memory was freed, @ref false if @c m was already cleared (or @c NULL)
bool mfile_free(mfile m);

/// @brief Reads the contents of the external file @c f, from the current position to the end of the file.
///
/// @param f the external file
///
/// @return a memory file with the contents of the external file, @c NULL if some error occurred
mfile mfile_read(FILE *f);

/// @brief Reads a segment of the external file @c f, of length @c len, starting at @c pos.
///
/// @param f the external file
/// @param pos the start position
/// @param len the length of the segment to read
///
/// @return a memory file with the contents of the external file, @c NULL if some error occurred
mfile mfile_read_from(FILE *f, ulong pos, ulong len);

/// @brief Reads the contents of the external file given by its path @c fileName.
///
/// @param fileName the path of the external file
///
/// @return a memory file with the contents of the external file, @c NULL if some error occurred
mfile mfile_load(char *fileName);

/// @brief Writes the entire content of the memory file @c m in the external file @c f, starting at the current position.
///
/// @param f the external file
/// @param m the memory file
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_write(FILE *f, mfile m);

/// @brief Writes a segment of the memory file @c m in the external file @c f, starting at position @c fpos.
///
/// @param f the external file
/// @param fpos the start position in the external file or @c -1 to start at the current position
/// @param m the memory file
/// @param mpos the start position in the memory file
/// @param len the length of the segment to write
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_write_to(FILE *f, long fpos, mfile m, ulong mpos, ulong len);

/// @brief Writes the entire content of the memory file @c m in the external file given by its path @c fileName.
///
/// @param fileName the path of the external file
/// @param m the memory file
/// @param append @ref true to append to the existing file, @ref false to replace its content
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_save(char *fileName, mfile m, bool append);

/// @brief Sets the position in the memory file @c m.
///
/// Same syntax as the system function @c fseek() for easy rewriting of the code.
///
/// @param m the memory file
/// @param pos the position
///
/// @return @ref true if successfull, @ref false otherwise
bool mseek(mfile m, ulong pos);

/// @brief Reads @c count blocks of size @c block from the memory file @c m to the given buffer @c buf.
///
/// Same syntax as the system function @c fread() for easy rewriting of the code.
///
/// @param buf the memory buffer to write to
/// @param block the size of one block
/// @param count the number of blocks to read
/// @param m the memory file
///
/// @return the number of blocks read
ulong mread(void *buf, ulong block, ulong count, mfile m);

/// @brief Writes @c count blocks of size @c block from the given buffer @c buf to the memory file @c m.
///
/// Same syntax as the system function @c fwrite() for easy rewriting of the code.
///
/// @param buf the memory buffer to read from
/// @param block the size of one block
/// @param count the number of blocks to write
/// @param m the memory file
///
/// @return the number of blocks written
ulong mwrite(void *buf, ulong block, ulong count, mfile m);

/// @brief Ensures that the memory buffer for the file @c m is at least @c the current position plus @c extra bytes.
///
/// @param m the memory file
/// @param extra the extra buffer capacity needed
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_ensure_cap(mfile m, ulong extra);

/// @brief Writes the byte @c c to the memory file @c m.
///
/// @param m the memory file
/// @param c the byte
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putb(mfile m, byte c);

/// @brief Writes the word @c w to the memory file @c m.
///
/// @param m the memory file
/// @param w the word
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putw(mfile m, word w);

/// @brief Writes the uint @c i to the memory file @c m.
///
/// @param m the memory file
/// @param i the uint
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_puti(mfile m, uint i);

/// @brief Writes the ulong @c l to the memory file @c m.
///
/// @param m the memory file
/// @param l the ulong
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putl(mfile m, ulong l);

/// @brief Writes the pointer @c p to the memory file @c m.
///
/// @param m the memory file
/// @param p the pointer
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putp(mfile m, void *p);

/// @brief Writes the double value @c d to the memory file @c m.
///
/// @param m the memory file
/// @param d the double value
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putd(mfile m, double d);

/// @brief Writes the @ref ldbl value @c l to the memory file @c m.
///
/// @param m the memory file
/// @param ld the @ref ldbl value
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putld(mfile m, ldbl ld);

/// @brief Writes a vector of bytes @c c to the memory file @c m.
///
/// @param m the memory file
/// @param c the vector of bytes
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putbs(mfile m, byte c[], ulong len);

/// @brief Writes a vector of words @c w to the memory file @c m.
///
/// @param m the memory file
/// @param w the vector of words
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putws(mfile m, word w[], ulong len);

/// @brief Writes a vector of uints @c i to the memory file @c m.
///
/// @param m the memory file
/// @param i the vector of uints
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putis(mfile m, uint i[], ulong len);

/// @brief Writes a vector of uints @c l to the memory file @c m.
///
/// @param m the memory file
/// @param l the vector of uints
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_putls(mfile m, ulong l[], ulong len);

/// @brief Reads a byte from the memory file @c m.
///
/// @param m the memory file
///
/// @return the byte, @c 0 if there was an error
byte mfile_getb(mfile m);

/// @brief Reads a word from the memory file @c m.
///
/// @param m the memory file
///
/// @return the word, @c 0 if there was an error
word mfile_getw(mfile m);

/// @brief Reads a uint from the memory file @c m.
///
/// @param m the memory file
///
/// @return the uint, @c 0 if there was an error
uint mfile_geti(mfile m);

/// @brief Reads a ulong from the memory file @c m.
///
/// @param m the memory file
///
/// @return the ulong, @c 0 if there was an error
ulong mfile_getl(mfile m);

/// @brief Reads a pointer from the memory file @c m.
///
/// @param m the memory file
///
/// @return the pointer, @c 0 if there was an error
void *mfile_getp(mfile m);

/// @brief Reads a double value from the memory file @c m.
///
/// @param m the memory file
///
/// @return the double value, @c NaN if there was an error
double mfile_getd(mfile m);

/// @brief Reads a @ref ldbl value from the memory file @c m.
///
/// @param m the memory file
///
/// @return the @ref ldbl value, @c NaN if there was an error
ldbl mfile_getld(mfile m);

/// @brief Reads a vector of bytes @c c from the memory file @c m.
///
/// @param m the memory file
/// @param c the vector of bytes
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_getbs(mfile m, byte c[], ulong len);

/// @brief Reads a vector of words @c w from the memory file @c m.
///
/// @param m the memory file
/// @param w the vector of words
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_getws(mfile m, word w[], ulong len);

/// @brief Reads a vector of uints @c i from the memory file @c m.
///
/// @param m the memory file
/// @param i the vector of uints
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_getis(mfile m, uint i[], ulong len);

/// @brief Reads a zero terminated string as a pointer in the memory file.
///
/// Increases the pointer in the memory file by the length of the string, including the null termination character.
/// Fails if there is no zero character starting from the current position.
///
/// @warning The returned string will be invalid if the memory file is freed !
///
/// @param m the memory file
///
/// @return the string at the current position, @ref NULL if some error occurred
char *mfile_get_str(mfile m);

/// @brief Reads a vector of ulongs @c l to the memory file @c m.
///
/// @param m the memory file
/// @param l the vector of ulongs
/// @param len the length of the vector
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_getls(mfile m, ulong l[], ulong len);

/// @brief Creates and returns a memory file containing the header of a binary file with given type ID, length of the
/// header and length of the entire file.
///
/// The position in the returned file is @c 20, after the @c 8 @c byte @c ID, the @c fileLength and the @c headerLength.
///
/// @param textId the ID of the type of binary file
/// @param fileLen the length of the entire file for which the header is created
/// @param headerLen the length of the header (and of the new memory file)
///
/// @return the memory file that contains the requested header, NULL if some error occurred
mfile mfile_header(char *textId, ulong fileLen, uint headerLen);

/// @brief Computes and adds the @c MD5 checksum to the end of the given @c header.
///
/// @param header the complete header, except the MD5 checksum
/// @param data the data in the file with the given header
/// @param dataLen the length of the @c data
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_header_add_md5(mfile header, void *data, ulong dataLen);

/// @brief Computes and checks the @c MD5 checksum of the given @c header.
///
/// @param header the complete header, containing the MD5 checksum as the last 16 bytes
/// @param data the data in the file with the given header
/// @param dataLen the length of the @c data
///
/// @return @ref true if the @c MD5 checksum could be verified, @ref false otherwise
bool mfile_header_check_md5(mfile header, void *data, ulong dataLen);

/// @brief Creates and returns a memory file containing the header of a binary file with given type ID, length of the
/// header and length of the entire file.
///
/// Similar to mfile_header(), but the reserved memory capacity is for the
/// entire file.
///
/// @param textId the ID of the type of binary file
/// @param fileLen the length of the entire file for which the header is created
/// @param headerLen the length of the header (and of the new memory file)
///
/// @return the memory file that contains the requested header, NULL if some error occurred
mfile mfile_from_header(char *textId, ulong fileLen, uint headerLen);

/// @brief Reads and returns a memory file containing the header of an external  binary file with given path
/// @c fileName and location.
///
/// It checks that the header is of a file with given type ID and its length is at least @c headerLen, fails if they do not match.
/// Leaves the pointer of the returned memory file at @c 0.
///
/// @param fileName the path of the binary file
/// @param pos the absolute position in the external file
/// @param textId the ID of the type of binary file
/// @param headerLen the minimum header lenght expected
///
/// @return the memory file that contains the requested header, @c NULL if some error occurred
mfile mfile_load_header(char *fileName, ulong pos, char *textId, int headerLen);

/// @brief Reads and returns a memory file containing the header of an external  binary file @c f at a given location.
///
/// @param f the binary file
/// @param pos the absolute position in the external file, negative to use the current position
///
/// @return the memory file that contains the requested header, @c NULL if some error occurred
mfile mfile_read_unknown_header(FILE *f, long pos);

/// @brief Checks if the @c header has the correct @c textId.
///
/// @param header the header
/// @param textId the text ID
///
/// @return @ref true if the header has the gien ID, @ref false otherwise
bool mfile_header_id(mfile header, char *textId);

/// @brief Reads and returns a memory file containing the header of an external binary file with given path
/// @c fileName.
///
/// @param fileName the path of the binary file
///
/// @return the memory file that contains the requested header, @c NULL if some error occurred
mfile mfile_load_unknown_header(char *fileName);

/// @brief Reads and returns a memory file containing the header of an external binary file @c f at a given location.
///
/// It checks that the header is of a file with given type ID and its length is at least @c headerLen, fails if they do not match.
/// Leaves the pointer of the file @c f after the end of the header and the pointer of the returned memory file at @c 0.
///
/// @param f the binary file
/// @param pos the absolute position in the external file, negative to use the current position
/// @param textId the ID of the type of binary file
/// @param headerLen the minimum header lenght expected
///
/// @return the memory file that contains the requested header, @c NULL if some error occurred
mfile mfile_read_header(FILE *f, long pos, char *textId, int headerLen);

/// @brief Creates, and writes to an existing memory file, a header of a binary file with given type ID, length of the
/// header and length of the entire file.
///
/// @param textId the ID of the type of binary file
/// @param fileLen the length of the entire file for which the header is created
/// @param headerLen the length of the header (and of the new memory file)
/// @param m the memory file to write to
///
/// @return @ref true if successfull, @ref false otherwise
bool mfile_write_header(char *textId, ulong fileLen, uint headerLen, mfile m);

/// @brief Checks that a correct header begins at position @c pos in the file of given type and name.
///
/// If so, it reads the entire file block described by the header in a @c mfile and returns it.
///
/// @param fileName the name of the file
/// @param pos the start position in the file
/// @param textId the type of the file
///
/// @return a memory file with the entire contents described by the header (may be shorter than the entire file)
mfile mfile_load_id(char *fileName, ulong pos, char *textId);

/// @brief Checks that a correct header begins at position @c pos in the given file of type @c textId.
///
/// If so, it reads the entire file block described by the header in a @c mfile and returns it.
///
/// @param f the the file
/// @param pos the start position in the file, negative to use the current position
/// @param textId the type of the file
///
/// @return a memory file with the entire contents described by the header (may be shorter than the entire file)
mfile mfile_read_id(FILE *f, long pos, char *textId);

#endif /* memFile_h */
