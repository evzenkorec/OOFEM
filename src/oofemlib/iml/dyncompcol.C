/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// inspired by SL++

#include "dyncompcol.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "mathfem.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"

namespace oofem {
REGISTER_SparseMtrx(DynCompCol, SMT_DynCompCol);

DynCompCol :: DynCompCol(void) : SparseMtrx(), base_(0)
{
#ifndef DynCompCol_USE_STL_SETS
    columns_ = NULL;
    rowind_  = NULL;
#else
    columns  = NULL;
#endif
}


DynCompCol :: DynCompCol(int n) : SparseMtrx(n, n), base_(0)
{
#ifndef DynCompCol_USE_STL_SETS
    columns_ = NULL;
    rowind_  = NULL;
#else
    columns  = NULL;
#endif

    nRows = nColumns = n;
}


/*****************************/
/*  Copy constructor         */
/*****************************/

DynCompCol :: DynCompCol(const DynCompCol &S) : SparseMtrx(S.nRows, S.nColumns), base_(S.base_)
{
#ifndef DynCompCol_USE_STL_SETS
    int i;
    if ( S.columns_ ) {
        this->columns_ = new FloatArray * [ S.nColumns ];
        for ( i = 0; i < S.nColumns; i++ ) {
            this->columns_ [ i ] = new FloatArray(*S.columns_ [ i ]);
        }
    } else {
        this->columns_ = NULL;
    }

    if ( S.rowind_ ) {
        this->rowind_ = new IntArray * [ S.nColumns ];
        for ( i = 0; i < S.nColumns; i++ ) {
            this->rowind_ [ i ] = new IntArray(*S.rowind_ [ i ]);
        }
    } else {
        this->rowind_ = NULL;
    }

#else
    int i;
    if ( S.columns ) {
        this->columns = new std :: map< int, double > * [ S.nColumns ];
        for ( i = 0; i < S.nColumns; i++ ) {
            this->columns [ i ] = new std :: map< int, double >(*S.columns [ i ]);
        }
    } else {
        this->columns = NULL;
    }

#endif

    this->nRows = S.nRows;
    this->nColumns = S.nColumns;
    this->version = S.version;
}


// Destructor
DynCompCol :: ~DynCompCol()
{
#ifndef DynCompCol_USE_STL_SETS
    int i;

    if ( this->columns_ ) {
        for ( i = 0; i < nColumns; i++ ) {
            delete this->columns_ [ i ];
        }

        delete this->columns_;
    }

    if ( this->rowind_ ) {
        for ( i = 0; i < nColumns; i++ ) {
            delete this->rowind_ [ i ];
        }

        delete this->rowind_;
    }

#else
    int i;

    if ( this->columns ) {
        for ( i = 0; i < nColumns; i++ ) {
            delete this->columns [ i ];
        }

        delete this->columns;
    }

#endif
}



/***************************/
/* Assignment operator...  */
/***************************/

DynCompCol &DynCompCol :: operator = ( const DynCompCol & C )
{
    base_   = C.base_;

#ifndef DynCompCol_USE_STL_SETS
    int i;

    if ( this->columns_ ) {
        for ( i = 0; i < nColumns; i++ ) {
            delete this->columns_ [ i ];
        }

        delete this->columns_;
    }

    if ( C.columns_ ) {
        this->columns_ = new FloatArray * [ C.nColumns ];
        for ( i = 0; i < C.nColumns; i++ ) {
            this->columns_ [ i ] = new FloatArray(*C.columns_ [ i ]);
        }
    } else {
        this->columns_ = NULL;
    }


    if ( this->rowind_ ) {
        for ( i = 0; i < nColumns; i++ ) {
            delete this->rowind_ [ i ];
        }

        delete this->rowind_;
    }

    if ( C.rowind_ ) {
        this->rowind_ = new IntArray * [ C.nColumns ];
        for ( i = 0; i < C.nColumns; i++ ) {
            this->rowind_ [ i ] = new IntArray(*C.rowind_ [ i ]);
        }
    } else {
        this->rowind_ = NULL;
    }

#else
    int i;

    if ( this->columns ) {
        for ( i = 0; i < nColumns; i++ ) {
            this->columns [ i ] = C.columns [ i ];
        }
    }

#endif

    nRows   = C.nRows;
    nColumns = C.nColumns;
    version = C.version;

    return * this;
}

SparseMtrx *DynCompCol :: GiveCopy() const
{
    DynCompCol *result = new DynCompCol(*this);
    return result;
}

void DynCompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    //      Check for compatible dimensions:
    if ( x.giveSize() != nColumns ) {
        OOFEM_ERROR("incompatible dimensions");
    }

    answer.resize(nRows);
    answer.zero();

#ifndef DynCompCol_USE_STL_SETS
    int j, t;
    double rhs;

    for ( j = 0; j < nColumns; j++ ) {
        rhs = x(j);
        for ( t = 1; t <= columns_ [ j ]->giveSize(); t++ ) {
            answer( rowind_ [ j ]->at(t) ) += columns_ [ j ]->at(t) * rhs;
        }
    }

#else
    double rhs;

    for ( int j = 0; j < nColumns; j++ ) {
        rhs = x(j);
        for ( auto &val: columns [ j ] ) {
            answer(val.first) += val.second * rhs;
        }
    }

#endif
}

void DynCompCol :: times(double x)
{
#ifndef DynCompCol_USE_STL_SETS
    for ( int j = 0; j < nColumns; j++ ) {
        columns_ [ j ]->times(x);
    }

#else
    for ( int j = 0; j < nColumns; j++ ) {
        for ( auto &val: columns [ j ] ) {
            val.second *= x;
        }
    }

#endif

    // increment version
    this->version++;
}

int DynCompCol :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    int neq = eModel->giveNumberOfDomainEquations(di, s);

#ifndef DynCompCol_USE_STL_SETS
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);

    nColumns = nRows = neq;

    if ( rowind_ ) {
        for ( int i = 0; i < nColumns; i++ ) {
            delete this->rowind_ [ i ];
        }

        delete this->rowind_;
    }

    rowind_ = ( IntArray ** ) new IntArray * [ neq ];
    for ( int j = 0; j < neq; j++ ) {
        rowind_ [ j ] = new IntArray();
    }

    // allocate value array
    if ( columns_ ) {
        for ( int i = 0; i < nColumns; i++ ) {
            delete this->columns_ [ i ];
        }

        delete this->columns_;
    }

    columns_ = ( FloatArray ** ) new FloatArray * [ neq ];
    for ( int j = 0; j < neq; j++ ) {
        columns_ [ j ] = new FloatArray();
    }

    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);

        for ( int ii : loc ) {
            if ( ii > 0 ) {
                for ( int jj : loc ) {
                    if ( jj > 0 ) {
                        this->insertRowInColumn(ii - 1, jj - 1);
                    }
                }
            }
        }
    }

    // loop over active boundary conditions
    std :: vector< IntArray >r_locs;
    std :: vector< IntArray >c_locs;

    for ( auto &gbc : domain->giveBcs() ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
        if ( bc != NULL ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                IntArray &krloc = r_locs [ k ];
                IntArray &kcloc = c_locs [ k ];
                for ( int ii : krloc ) {
                    if ( ii > 0 ) {
                        for ( int jj : kcloc ) {
                            if ( jj > 0 ) {
                                this->insertRowInColumn(jj - 1, ii - 1);
                            }
                        }
                    }
                }
            }
        }
    }

    int nz_ = 0;
    for ( int j = 0; j < neq; j++ ) {
        nz_ += this->rowind_ [ j ]->giveSize();
    }

    OOFEM_LOG_DEBUG("DynCompCol info: neq is %d, nelem is %d\n", neq, nz_);
#else
    nColumns = nRows = neq;

    columns.resize( neq );
    for ( auto &col: columns ) {
        col.clear();
    }

#endif

    // increment version
    this->version++;

    return true;
}


int DynCompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;

#  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'k' and 'loc' mismatch");
    }

#  endif

    dim = mat.giveNumberOfRows();

    for ( j = 1; j <= dim; j++ ) {
        jj = loc.at(j);
        if ( jj ) {
            for ( i = 1; i <= dim; i++ ) {
                ii = loc.at(i);
                if ( ii ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

int DynCompCol :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
#ifndef DynCompCol_USE_STL_SETS
    /*
     * /// slow and safe
     * int        i,j,ii,jj,dim1,dim2 ;
     *
     * this->checkSizeTowards(rloc, cloc);
     *
     * dim1 = mat.giveNumberOfRows() ;
     * dim2 = mat.giveNumberOfColumns() ;
     * for (i=1 ; i<= dim1; i++) {
     * ii = rloc.at(i);
     * if (ii) for (j=1 ; j<= dim2; j++) {
     * jj = cloc.at(j);
     * if (jj) this->at(ii,jj) += mat.at(i,j);
     * }
     * }
     * return 1;
     */

    /// optimized low-end implementation
    IntArray rowsToAdd( rloc.giveSize() );
    int i, ii, ii1, j, jj, jj1, rowindx;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    /*
     * // adjust the size of receiver
     * int maxid=0;
     * for (i=0; i<rsize; i++) maxid = max(maxid, rloc(i));
     * for (i=0; i<csize; i++) maxid = max(maxid, cloc(i));
     * this->growTo (maxid);
     */

    for ( i = 0; i < csize; i++ ) {
        if ( ( ii = cloc(i) ) ) {
            ii1 = ii - 1;
            for ( j = 0; j < rsize; j++ ) {
                if ( ( jj = rloc(j) ) ) {
                    jj1 = jj - 1;

                    /*
                     *       rowindx = 0;
                     *   for (int t=columns_[ii1]->giveSize(); t > 0; t--)
                     *    if (rowind_[ii1]->at(t) == (jj1)) {rowindx = t; break;}
                     */

                    //     rowindx = this->giveRowIndx (ii1, jj1);
                    //     if (rowindx == 0) {
                    /*
                     *    int oldsize = rowind_[ii1]->giveSize();
                     *    rowind_[ii1]->resizeWithValues(oldsize+1, DynCompCol_CHUNK);
                     *    columns_[ii1]->resizeWithValues(oldsize+1, DynCompCol_CHUNK);
                     *    rowindx = oldsize+1;
                     *    rowind_[ii1]->at(oldsize+1) = jj1;
                     */
                    //      rowindx = this->insertRowInColumn (ii1, jj1);
                    //     }

                    rowindx = this->insertRowInColumn(ii1, jj1);
                    //if (rowind_[j-1]->at(t) == (i-1)) return columns_[j-1]->at(t);
                    columns_ [ ii1 ]->at(rowindx) += mat(j, i);
                }
            }
        }
    }

#else
    int i, ii, ii1, j, jj, jj1;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    for ( i = 0; i < csize; i++ ) {
        if ( ( ii = cloc(i) ) ) {
            ii1 = ii - 1;
            for ( j = 0; j < rsize; j++ ) {
                if ( ( jj = rloc(j) ) ) {
                    jj1 = jj - 1;
                    this->columns [ ii1 ] [ jj1 ] += mat(j, i);
                }
            }
        }
    }

#endif
    // increment version
    this->version++;

    return 1;
}

void DynCompCol :: zero()
{
#ifndef DynCompCol_USE_STL_SETS
    for ( int j = 0; j < nColumns; j++ ) {
        columns_ [ j ]->zero();
    }

#else
    for ( auto &col: columns ) {
        for ( auto &val: col ) {
            val.second = 0.;
        }
    }

#endif
    // increment version
    this->version++;
}

void DynCompCol :: printStatistics() const
{
    int j, nz_ = 0;
#ifndef DynCompCol_USE_STL_SETS
    for ( j = 0; j < nColumns; j++ ) {
        nz_ += rowind_ [ j ]->giveSize();
    }

#else
    for ( j = 0; j < nColumns; j++ ) {
        nz_ += columns [ j ].size();
#endif
    OOFEM_LOG_DEBUG("DynCompCol info: neq is %d, nelem is %d\n", nColumns, nz_);
}


/*********************/
/*   Array access    */
/*********************/

double &DynCompCol :: at(int i, int j)
{
    // increment version
    this->version++;
#ifndef DynCompCol_USE_STL_SETS
    /*
     * for (int t=1; t<=columns_[j-1]->giveSize(); t++)
     * if (rowind_[j-1]->at(t) == (i-1)) return columns_[j-1]->at(t);
     * printf ("DynCompCol::operator(): Array accessing exception -- out of bounds.\n");
     * exit(1);
     * return columns_[0]->at(1);   // return to suppress compiler warning message
     */
    int rowIndx;
    if ( ( rowIndx = this->giveRowIndx(j - 1, i - 1) ) ) {
        return columns_ [ j - 1 ]->at(rowIndx);
    }

    OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
    return columns_ [ 0 ]->at(1); // return to suppress compiler warning message

#else
    return this->columns [ j - 1 ] [ i - 1 ];

#endif
}


double DynCompCol :: at(int i, int j) const
{
#ifndef DynCompCol_USE_STL_SETS
    /*
     * for (int t=1; t<=columns_[j-1]->giveSize(); t++)
     * if (rowind_[j-1]->at(t) == (i-1)) return columns_[j-1]->at(t);
     * if (i <= nRows && j <= nColumns) return 0.0;
     * else
     * {
     * printf ("DynCompCol::operator(): Array accessing exception -- out of bounds.\n");
     * exit(1);
     * return (0);   // return to suppress compiler warning message
     * }
     */
    int rowIndx;
    if ( ( rowIndx = this->giveRowIndx(j - 1, i - 1) ) ) {
        return columns_ [ j - 1 ]->at(rowIndx);
    }

    if ( i <= nRows && j <= nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return columns_ [ 0 ]->at(1); // return to suppress compiler warning message
    }

#else
    auto pos = this->columns [ j - 1 ].find(i - 1);
    if ( pos != this->columns [ j - 1 ].end() ) {
        return pos->second;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return 0.0;
    }

#endif
}

double DynCompCol :: operator() (int i, int j)  const
{
#ifndef DynCompCol_USE_STL_SETS
    /*
     * for (int t=1; t<=columns_[j]->giveSize(); t++)
     * if (rowind_[j]->at(t) == i) return columns_[j]->at(t);
     * if (i < nRows && j < nColumns) return 0.0;
     * else
     * {
     * printf ("DynCompCol::operator(): Array accessing exception -- out of bounds.\n");
     * exit(1);
     * return (0);   // return to suppress compiler warning message
     * }
     */
    int rowIndx;
    if ( ( rowIndx = this->giveRowIndx(j, i) ) ) {
        return columns_ [ j ]->at(rowIndx);
    }

    if ( i < nRows && j < nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return columns_ [ 0 ]->at(1); // return to suppress compiler warning message
    }

#else
    auto pos = this->columns [ j ].find(i);
    if ( pos != this->columns [ j ].end() ) {
        return pos->second;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return 0.0;
    }

#endif
}

double &DynCompCol :: operator() (int i, int j)
{
    // increment version
    this->version++;
#ifndef DynCompCol_USE_STL_SETS
    /*
     * for (int t=1; t<=columns_[j]->giveSize(); t++)
     * if (rowind_[j]->at(t) == i) return columns_[j]->at(t);
     *
     * printf ("DynCompCol::set: Array element (%d,%d) not in sparse structure -- cannot assign.\n",i,j);
     * exit(1);
     * return columns_[j]->at(1);   // return to suppress compiler warning message
     */

    int rowIndx;
    if ( ( rowIndx = this->giveRowIndx(j, i) ) ) {
        return columns_ [ j ]->at(rowIndx);
    }

    OOFEM_ERROR("Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return columns_ [ 0 ]->at(1); // return to suppress compiler warning message

#else
    return this->columns [ j ] [ i ];

#endif
}


void DynCompCol :: timesT(const FloatArray &x, FloatArray &answer) const
{
    // Check for compatible dimensions:
    if ( x.giveSize() != nRows ) {
        OOFEM_ERROR("Error in CompCol -- incompatible dimensions");
    }
    answer.resize(nColumns);
    answer.zero();

#ifndef DynCompCol_USE_STL_SETS
    int i, t;
    double r;

    for ( i = 0; i < nColumns; i++ ) {
        r = 0.0;
        for ( t = 1; t <= columns_ [ i ]->giveSize(); t++ ) {
            r += columns_ [ i ]->at(t) * x( rowind_ [ i ]->at(t) );
        }

        answer(i) = r;
    }

#else
    double r;
    for ( int i = 0; i < nColumns; i++ ) {
        double r = 0.0;
        for ( auto &val: columns [ i ] ) {
            r += val.second * x(val.first);
        }

        answer(i) = r;
    }

#endif
}



void DynCompCol :: checkSizeTowards(IntArray &loc)
{
    int i, maxid = 0;
    int size = loc.giveSize();
    // adjust the size of receiver
    for ( i = 0; i < size; i++ ) {
        maxid = max( maxid, loc(i) );
    }

    this->growTo(maxid);

#ifndef DynCompCol_USE_STL_SETS
    int ii, j, jj;

    for ( i = 0; i < size; i++ ) {
        if ( ( ii = loc(i) ) ) {
            for ( j = 0; j < size; j++ ) {
                if ( ( jj = loc(j) ) ) {
                    this->insertRowInColumn(ii - 1, jj - 1);
                }
            }
        }
    }

#endif
}


void DynCompCol :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
{
    int i, maxid = 0;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    // adjust the size of receiver
    for ( i = 0; i < rsize; i++ ) {
        maxid = max( maxid, rloc(i) );
    }

    for ( i = 0; i < csize; i++ ) {
        maxid = max( maxid, cloc(i) );
    }

    this->growTo(maxid);

#ifndef DynCompCol_USE_STL_SETS
    IntArray rowsToAdd( rloc.giveSize() );
    int ii, j, jj;

    for ( i = 0; i < csize; i++ ) {
        if ( ( ii = cloc(i) ) ) {
            for ( j = 0; j < rsize; j++ ) {
                if ( ( jj = rloc(j) ) ) {
                    insertRowInColumn(ii - 1, jj - 1);
                }
            }
        }
    }

#endif
}



void DynCompCol :: growTo(int ns)
{
    if ( ns > nColumns ) {
#ifndef DynCompCol_USE_STL_SETS
        FloatArray **newcolumns_ = new FloatArray * [ ns ];
        IntArray **newrowind_ = new IntArray * [ ns ];

        // copy existing columns
        for ( int i = 0; i < nColumns; i++ ) {
            newcolumns_ [ i ] = columns_ [ i ];
            newrowind_ [ i ] = rowind_ [ i ];
        }

        delete columns_;
        delete rowind_;

        columns_ = newcolumns_;
        rowind_  = newrowind_;
#else
        columns.resize(ns);
#endif

        nColumns = nRows = ns;
    }
}


#ifndef DynCompCol_USE_STL_SETS
int DynCompCol :: giveRowIndx(int col, int row) const
{
    // fast row indx search, based on assumption, that row indices are sorted
    int left = 1, right = this->rowind_ [ col ]->giveSize();
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( right == 0 ) {
        return 0;
    }

    if ( this->rowind_ [ col ]->at(right) == row ) {
        return right;
    }

    while ( !( ( ( middleVal = this->rowind_ [ col ]->at(middle) ) == row ) || ( middle == left ) ) ) {
        if ( row > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == row ) {
        return middle;
    }

    return 0;
}



int
DynCompCol :: insertRowInColumn(int col, int row)
{
    // insert row into column, preserving order of row indexes.
    int i, oldsize = this->rowind_ [ col ]->giveSize();
    int left = 1, right = oldsize;
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( oldsize == 0 ) {
        rowind_ [ col ]->resizeWithValues(1, DynCompCol_CHUNK);
        columns_ [ col ]->resizeWithValues(1, DynCompCol_CHUNK);
        columns_ [ col ]->at(1) = 0.0;
        rowind_ [ col ]->at(1) = row;
        return 1;
    }

    if ( this->rowind_ [ col ]->at(right) == row ) {
        return right;
    }

    while ( !( ( ( middleVal = this->rowind_ [ col ]->at(middle) ) == row ) || ( middle == left ) ) ) {
        if ( row > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == row ) {
        return middle;
    }

    // we have to insert new row entry
    if ( row > this->rowind_ [ col ]->at(oldsize) ) {
        right = oldsize + 1;
    } else if ( row < this->rowind_ [ col ]->at(1) ) {
        right = 1;
    }

    // insert row at middle+1 position
    rowind_ [ col ]->resizeWithValues(oldsize + 1, DynCompCol_CHUNK);
    columns_ [ col ]->resizeWithValues(oldsize + 1, DynCompCol_CHUNK);

    for ( i = oldsize; i >= right; i-- ) {
        rowind_ [ col ]->at(i + 1) = rowind_ [ col ]->at(i);
        columns_ [ col ]->at(i + 1) = columns_ [ col ]->at(i);
    }

    columns_ [ col ]->at(right) = 0.0;
    rowind_ [ col ]->at(right) = row;
    return right;
}

#endif
} // end namespace oofem
