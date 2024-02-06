//
//  treeMap.h
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
 \file treeMap.h
 \brief A data structure and a collection of basic functions to manipulate dyadic trees.
 
 A tree is the representation of a set or of a positive measure in the plane, with a limited precision.
 If the original set contains points that are closer together than the resolution of the tree, only their
 count will be retained in the corresponding dyadic square.
 
 @see planarSet.h
 @see nSet.h
 @see dRect.h
 @see bitmap.h 
*/

#ifndef treeMap_h
#define treeMap_h

#include <stdio.h>
#include <stdlib.h>

#include "nSet.h"
#include "fp80.h"
#include "memFile.h"
#include "dRect.h"
#include "bitmap.h"

/// Max level step from a parent node to child node.
#define TMAP_MAX_LEVEL_STEP      6

/// Root level to represent the Mandelbrot set.
#define TMAP_MANDEL_MIN_LEVEL   -2

/// X coord of the root node to represent the Mandelbrot set.
#define TMAP_MANDEL_X            0

/// Y coord of the root node to represent the Mandelbrot set.
#define TMAP_MANDEL_Y            0

/// Max level of a leaf.
#define TMAP_MAX_LEVEL          30

/// Binary file ID for a serialized \ref treeMap.
#define TMAP_FILE_ID  "tmap v01"

/// Binary file header length of a serialized @ref tmap.
#define TMAP_HEADER_LEN  68

struct treeMap_struct;

/// A pointer to a treeMap_struct.
typedef struct treeMap_struct *tmap;

/// @union mapCount_union
/// \brief A pointer to a child node or the count of points present in a sub-square of a leaf.
///
/// This is for efficient use of memory, while having a single data structure for both types of nodes.
typedef union {
    tmap ptr;
    ulong count;
} mapCount_union;

/// @struct treeMap_struct
/// \brief This data structure caracterizes and stores a dyadic tree.
///
/// Leafs contain the number of points that are present in each sub-square, using the same memory
/// that is reserved for pointers to children (when the node is not a leaf). Nodes are not aware
/// of their parents, therefore operations should be performed only on the root.
typedef struct treeMap_struct {
    ulong count;      ///< he number of points present in the represented square
    char level;       ///< the level of the node, the size of the square is @c 2^{-level}
    char maxLevel;    ///< max level of any leaf which has this node as an ancestor
    char levelStep;   ///< the level difference between this node and its children
    char minLevel;    ///< the level of the root
    word children;    ///< the number of [non-empty] children of this node
    int x;            ///< absolute coordinate of the left side of the square, to be divided by @c 2^{level}
    int y;            ///< absolute coordinate of the bottom side of the square, to be divided by @c 2^{level}
    mapCount_union chld[]; ///< the list of children if this node is not a leaf, the count of points in each sub-square otherwise
} treeMap_struct;

/// Tests if @c t is a root.
#define tmap_isRoot(t) (t->level == t->minLevel)

/// Tests if @c t is a leaf.
#define tmap_isLeaf(t) (t->level + t->levelStep > t->maxLevel)

/// The maximal number of children of the node @c t.
#define tmap_maxChildren(t) (1 << (t->levelStep << 1))

/// The number of columns of the node @c t.
#define tmap_columns(t) (1 << t->levelStep)

/// The number bytes need to store the positions of the children of @c t.
#define tmap_bitDataLen(t) ((tmap_maxChildren(t) - 1) / 8 + 1)

/// The maximal number bytes need to store all the counts of the leaf @c t.
#define tmap_maxLeafDataLen(t) (8 * tmap_maxChildren(t))

/// Max (exclusive) of number of points in a sub-square of a leaf that are to be stored in one byte in th binary file.
#define TMAP_MAX_CHAR    (1L <<  7)
/// Max (exclusive) of number of points in a sub-square of a leaf that are to be stored in two bytes in th binary file.
#define TMAP_MAX_SHORT   (1L << 14)
/// Max (exclusive) of number of points in a sub-square of a leaf that are to be stored in four bytes in th binary file.
#define TMAP_MAX_INT     (1L << 29)

// MARK: creation and deletion of nodes

/// @brief Creates and returns a new node of a tree.
///
/// @param x the @c x coordiante of the left of the square represented by this node
/// @param y the @c y coordiante of the bottom of the square represented by this node
/// @param level the level of the node, the size of the represented square is @c 2^{-level}
/// @param maxLevel the max level of a node in this tree
/// @param levelStep the difference of levels between a parent node and a child node
/// @param minLevel the level of the root
///
/// @return the new node, @c NULL if there is an error in the parameters
tmap tmap_new(int x, int y, char level, char maxLevel, char levelStep, char minLevel);

/// @brief Creates and returns a new root of a tree representing points in the Mandelbrot set.
///
/// @param maxLevel the max level of a node in this tree
/// @param levelStep the difference of levels between a parent node and a child node
///
/// @return the new node, @c NULL if there is an error in the parameters
tmap tmap_new_mandel(char maxLevel, char levelStep);

/// @brief Computes and returns the map of the given \ref nset of points included in the Mandelbrot set.
///
/// @param ps the set of points
/// @param maxLevel the max level of nodes in the tree
/// @param levelStep the level difference of parent and children nodes in the tree
///
/// @return the new tree, @c NULL if there is an error in the parameters
tmap tmap_map(nset_t ps, char maxLevel, char levelStep);

/// @brief Frees the memory used by this node and all its children.
///
/// @param tm the node
void tmap_free(tmap tm);

/// @brief Creates and return a copy of the tree which has this node as a root.
///
/// @param tree the node
///
/// @return the new tree, @c NULL is @c tree is @c NULL
tmap tmap_clone(tmap tree);

/// @brief Creates and adds a new child node to @c tm, at the position @c (i,j).
///
/// Fails if there is already a child on this position.
///
/// @param tm the node
/// @param i the column of @c tm, starting at @c 0
/// @param j the row of @c tm, starting at @c 0
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_new_child(tmap tm, int i, int j);

// MARK: addition of points, sets of points and trees

/// @brief Adds a point to the tree.
///
/// @param tree the [root of the] tree
/// @param p the point
///
/// @return @ref true if the point was added,  @ref false otherwise
bool tmap_add(tmap tree, u128 p);

/// @brief Adds a point to the tree.
///
/// @param tree the [root of the] tree
/// @param x the upper 64-bit of the \ref uint128 @c x coordinate of the point
/// @param y the upper 64-bit of the \ref uint128 @c y coordinate of the point
///
/// @return @ref true if the point was added,  @ref false otherwise
bool tmap_addl(tmap tree, ulong x, ulong y);

/// @brief Adds a set of points @c ps to the @c tree.
///
/// @param tree the tree
/// @param ps the set of points
///
/// @return the number of points that were added
ulong tmap_adds(tmap tree, nset_t ps);

/// @brief Computes the union of the two trees and stores the result in @c tree.
///
/// @param tree the destination tree
/// @param op some other
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_union(tmap tree, tmap op);

/// @brief Reduces the max level of the given @c tree.
///
/// @param tree the tree
/// @param maxLevel the new maxLevel
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_reduce(tmap tree, int maxLevel);

// MARK: comparisons and node counts

/// @brief Compares two trees.
///
/// @param t1 a tree
/// @param t2 another tree
///
/// @return @ref true if the trees represent the same set / measure, @ref false otherwise
bool tmap_eq(tmap t1, tmap t2);

/// @brief Compares two trees.
///
/// @param t1 a tree
/// @param t2 another tree
///
/// @return @ref true if @c t1<=t2, @ref false otherwise
bool tmap_leq(tmap t1, tmap t2);

/// @brief Counts the number of nodes of the tree @c tm.
///
/// @param tm the root of the tree.
///
/// @return the number of nodes
ulong tmap_nodes(tmap tm);

/// @brief Return the memory used by the tree @c tm.
///
/// @param tm the root of the tree.
///
/// @return the number of bytes of memory used
ulong tmap_memory(tmap tm);

/// @brief Returns the number of points of the tree that lay in the given dyadic rectangle relative coordinates.
///
/// @param tm the tree
/// @param r the dyadic rectangle
///
/// @return the number of points of the tree that lay in the given dyadic rectangle
ulong tmap_count(tmap tm, drect r);

/// @brief Returns the number of points of the tree that lay in the given dyadic rectangle with absolute coordinates.
///
/// @param tm the tree
/// @param r the dyadic rectangle
///
/// @return the number of points of the tree that lay in the given dyadic rectangle
ulong tmap_count_abs(tmap tm, drect r);

// MARK: bitmap operations

/// @brief Adds the count of points in the tree for each dyadic square represented by a pixel of the given bitmap @c bm.
///
/// @param bm the bitmap
/// @param tm the tree
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_add_counts(bmap bm, tmap tm);

/// @brief Computes the smalles dyadic rectangle that contains the given tree @c tm.
///
/// @param boundingRect the result
/// @param tm the tree
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_bounds(drect boundingRect, tmap tm);

/// @brief Return a bitmap of the tree in the given dyadic rectangle @c r.
///
/// @param tm the tree
/// @param r the dyadic rectangle
///
/// @return the bitmap, @c NULL if some error occurred
bmap tmap_bitmap(tmap tm, drect r);

// MARK: IO operations

/// @brief Reads and returns a tree from the beginning of the memory file @c m.
///
/// @param m the @ref mfile name
/// @param pos the absolute position in the file, negative to use the current position
/// @param checkMD5 @ref true to check the MD5 sum of the file
/// @param bounds the bounds of the tree stored in the file, unused if @c NULL
///
/// @return the tree stored in the file, @c NULL if some error occurred
tmap tmap_mread(mfile m, bool checkMD5, drect bounds);

/// @brief Reads and returns a tree from the file @c f.
///
/// @param f the file
/// @param pos the absolute position in the file, negative to use the current position
/// @param checkMD5 @ref true to check the MD5 sum of the file
/// @param bounds the bounds of the tree stored in the file, unused if @c NULL
///
/// @return the tree stored in the file, @c NULL if some error occurred
tmap tmap_read(FILE *f, long pos, bool checkMD5, drect bounds);

/// @brief Reads and returns a tree from the file with given @c fileName.
///
/// @param fileName the file name
/// @param pos the absolute position in the file
/// @param checkMD5 @ref true to check the MD5 sum of the file
/// @param bounds the bounds of the tree stored in the file, unused if @c NULL
///
/// @return the tree stored in the file, @c NULL if some error occurred
tmap tmap_load(char *fileName, ulong pos, bool checkMD5, drect bounds);

/// @brief Writes the tree to a new memory file.
///
/// Fails if the node @c tree is not the root.
///
/// @param tree the tree
/// @param bounds the bounds of the tree, @c NULL to compute on the fly
///
/// @return the memory file, @c NULL if some error occurred
mfile tmap_write(tmap tree, drect bounds);

/// @brief Writes the tree to the file with path @c fileName.
///
/// Fails if the node @c tree is not the root.
///
/// @param tree the tree
/// @param fileName the file name
///
/// @return @ref true if successfull, @ref false otherwise
bool tmap_save(tmap tree, char *fileName);

#endif /* treeMap_h */
