/*****************************************************************************/
/* Name: matrix.h                                                            */
/* Uses: Class for matrix math functions.                                    */
/* Date: 4/19/2011                                                           */
/* Author: Andrew Que <http://www.DrQue.net/>                                */
/* Revisions:                                                                */
/*   0.1 - 2011/04/19 - QUE - Creation.                                      */
/*   0.5 - 2011/04/24 - QUE - Most functions are complete.                   */
/*   0.8 - 2011/05/01 - QUE -                                                */
/*     = Bug fixes.                                                          */
/*     + Dot product.                                                        */
/*   1.0 - 2011/11/26 - QUE - Release.                                       */
/*                                                                           */
/* Notes:                                                                    */
/*   This unit implements some very basic matrix functions, which include:   */
/*    + Addition/subtraction                                                 */
/*    + Transpose                                                            */
/*    + Row echelon reduction                                                */
/*    + Determinant                                                          */
/*    + Dot product                                                          */
/*    + Matrix product                                                       */
/*    + Scalar product                                                       */
/*    + Inversion                                                            */
/*    + LU factorization/decomposition                                       */
/*     There isn't much for optimization in this unit as it was designed as  */
/*   more of a learning experience.                                          */
/*                                                                           */
/* License:                                                                  */
/*   This program is free software: you can redistribute it and/or modify    */
/*   it under the terms of the GNU General Public License as published by    */
/*   the Free Software Foundation, either version 3 of the License, or       */
/*   (at your option) any later version.                                     */
/*                                                                           */
/*   This program is distributed in the hope that it will be useful,         */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*   GNU General Public License for more details.                            */
/*                                                                           */
/*   You should have received a copy of the GNU General Public License       */
/*   along with this program.  If not, see <http://www.gnu.org/licenses/>.   */
/*                                                                           */
/*                     (C) Copyright 2011 by Andrew Que                      */
/*                           http://www.DrQue.net/                           */
/*                                                                           */
/*****************************************************************************/
/* Note (P.F. Marteau): some slight modifications have been added to comply  */
/* with the matrix template type defined in the dlib library c.f.            */
/*                                            www.http://dlib.net/           */
/*****************************************************************************/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <fstream>
#include <cassert>
#include <climits>
#include <vector>

// Class forward for identity matrix.
template< class TYPE > class IdentityMatrix;

//=============================================================================
// Matrix template class
//   Contains a set of matrix manipulation functions.  The template is designed
// so that the values of the matrix can be of any type that allows basic
// arithmetic.
//=============================================================================
template< class TYPE = int >
  class matrix
  {
    protected:
      // Matrix data.
      unsigned rows;
      unsigned columns;

      // Storage for matrix data.
      std::vector< std::vector< TYPE > > mat;

      // Order sub-index for rows.
      //   Use: matrix[ order[ row ] ][ column ].
      unsigned * order;

      //-------------------------------------------------------------
      // Return the number of leading zeros in the given row.
      //-------------------------------------------------------------
      unsigned getLeadingZeros
      (
        // Row to count
        unsigned row
      ) const
      {
        TYPE const ZERO = static_cast< TYPE >( 0 );
        unsigned column = 0;
        while ( ZERO == mat[ row ][ column ] )
          ++column;

        return column;
      }

      //-------------------------------------------------------------
      // Reorder the matrix so the rows with the most zeros are at
      // the end, and those with the least at the beginning.
      //
      // NOTE: The matrix data itself is not manipulated, just the
      // 'order' sub-indexes.
      //-------------------------------------------------------------
      void reorder()
      {
        unsigned * zeros = new unsigned[ rows ];

        for ( unsigned row = 0; row < rows; ++row )
        {
          order[ row ] = row;
          zeros[ row ] = getLeadingZeros( row );
        }

        for ( unsigned row = 0; row < (rows-1); ++row )
        {
          unsigned swapRow = row;
          for ( unsigned subRow = row + 1; subRow < rows; ++subRow )
          {
            if ( zeros[ order[ subRow ] ] < zeros[ order[ swapRow ] ] )
              swapRow = subRow;
          }

          unsigned hold    = order[ row ];
          order[ row ]     = order[ swapRow ];
          order[ swapRow ] = hold;
        }

        delete zeros;
      }

      //-------------------------------------------------------------
      // Divide a row by given value.  An elementary row operation.
      //-------------------------------------------------------------
      void divideRow
      (
        // Row to divide.
        unsigned row,

        // Divisor.
        TYPE const & divisor
      )
      {
        for ( unsigned column = 0; column < columns; ++column )
          mat[ row ][ column ] /= divisor;
      }

      //-------------------------------------------------------------
      // Modify a row by adding a scaled row. An elementary row
      // operation.
      //-------------------------------------------------------------
      void rowOperation
      (
        unsigned row,
        unsigned addRow,
        TYPE const & scale
      )
      {
        for ( unsigned column = 0; column < columns; ++column )
          mat[ row ][ column ] += mat[ addRow ][ column ] * scale;
      }

      //-------------------------------------------------------------
      // Allocate memory for matrix data.
      //-------------------------------------------------------------
      void allocate
      (
        unsigned rowNumber,
        unsigned columnNumber
      )
      {
        // Allocate order integers.
        order = new unsigned[ rowNumber ];

        // Setup matrix sizes.
        mat.resize( rowNumber );
        for ( unsigned row = 0; row < rowNumber; ++row )
          mat[ row ].resize( columnNumber );
      }

      //-------------------------------------------------------------
      // Free memory used for matrix data.
      //-------------------------------------------------------------
      void deallocate
      (
        unsigned rowNumber,
        unsigned columnNumber
      )
      {
        // Free memory used for storing order (if there is any).
        if ( 0 != rowNumber )
          delete[] order;
      }

    public:
      // Used for matrix concatenation.
      typedef enum
      {
        TO_RIGHT,
        TO_BOTTOM
      } Position;

      //-------------------------------------------------------------
      // Return the number of rows in this matrix.
      //-------------------------------------------------------------
      unsigned getRows() const
      {
        return rows;
      }
      unsigned nr() const
      {
        return rows;
      }

      //-------------------------------------------------------------
      // Return the number of columns in this matrix.
      //-------------------------------------------------------------
      unsigned getColumns() const
      {
        return columns;
      }
      unsigned nc() const
      {
        return columns;
      }

      //-------------------------------------------------------------
      // Get an element of the matrix.
      //-------------------------------------------------------------
      TYPE get
      (
        unsigned row,   // Which row.
        unsigned column // Which column.
      ) const
      {
        assert( row < rows );
        assert( column < columns );

        return mat[ row ][ column ];
      }

      //-------------------------------------------------------------
      // Proform LU decomposition.
      // This will create matrices L and U such that A=LxU
      //-------------------------------------------------------------
      void LU_Decomposition
      (
        matrix & upper,
        matrix & lower
      ) const
      {
        assert( rows == columns );

        TYPE const ZERO = static_cast< TYPE >( 0 );

        upper = *this;
        lower = *this;

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            lower.mat[ row ][ column ] = ZERO;

        for ( unsigned row = 0; row < rows; ++row )
        {
          TYPE value = upper.mat[ row ][ row ];
          if ( ZERO != value )
          {
            upper.divideRow( row, value );
            lower.mat[ row ][ row ] = value;
          }

          for ( unsigned subRow = row + 1; subRow < rows; ++subRow )
          {
            TYPE value = upper.mat[ subRow ][ row ];
            upper.rowOperation( subRow, row, -value );
            lower.mat[ subRow ][ row ] = value;
          }
        }
      }

      //-------------------------------------------------------------
      // Set an element in the matrix.
      //-------------------------------------------------------------
      void put
      (
        unsigned row,
        unsigned column,
        TYPE const & value
      )
      {
        assert( row < rows );
        assert( column < columns );

        mat[ row ][ column ] = value;
      }

      //-------------------------------------------------------------
      // Return part of the matrix.
      // NOTE: The end points are the last elements copied.  They can
      // be equal to the first element when wanting just a single row
      // or column.  However, the span of the total matrix is
      // ( 0, rows - 1, 0, columns - 1 ).
      //-------------------------------------------------------------
      matrix getSubMatrix
      (
        unsigned startRow,
        unsigned endRow,
        unsigned startColumn,
        unsigned endColumn,
        unsigned const * newOrder = NULL
      )
      {
        matrix subMatrix( endRow - startRow + 1, endColumn - startColumn + 1 );

        for ( unsigned row = startRow; row <= endRow; ++row )
        {
          unsigned subRow;
          if ( NULL == newOrder )
            subRow = row;
          else
            subRow = newOrder[ row ];

          for ( unsigned column = startColumn; column <= endColumn; ++column )
            subMatrix.mat[ row - startRow ][ column - startColumn ] =
              mat[ subRow ][ column ];
        }

        return subMatrix;
      }

      //-------------------------------------------------------------
      // Return a single column from the matrix.
      //-------------------------------------------------------------
      matrix getColumn
      (
        unsigned column
      )
      {
        return getSubMatrix( 0, rows - 1, column, column );
      }

      //-------------------------------------------------------------
      // Return a single row from the matrix.
      //-------------------------------------------------------------
      matrix getRow
      (
        unsigned row
      )
      {
        return getSubMatrix( row, row, 0, columns - 1 );
      }

      //-------------------------------------------------------------
      // Place matrix in reduced row echelon form.
      //-------------------------------------------------------------
      void reducedRowEcholon()
      {
        TYPE const ZERO = static_cast< TYPE >( 0 );

        // For each row...
        for ( unsigned rowIndex = 0; rowIndex < rows; ++rowIndex )
        {
          // Reorder the rows.
          reorder();

          unsigned row = order[ rowIndex ];

          // Divide row down so first term is 1.
          unsigned column = getLeadingZeros( row );
          TYPE divisor = mat[ row ][ column ];
          if ( ZERO != divisor )
          {
            divideRow( row, divisor );

            // Subtract this row from all subsequent rows.
            for ( unsigned subRowIndex = ( rowIndex + 1 ); subRowIndex < rows; ++subRowIndex )
            {
              unsigned subRow = order[ subRowIndex ];
              if ( ZERO != mat[ subRow ][ column ] )
                rowOperation
                (
                  subRow,
                  row,
                  -mat[ subRow ][ column ]
                );
            }
          }

        }

        // Back substitute all lower rows.
        for ( unsigned rowIndex = ( rows - 1 ); rowIndex > 0; --rowIndex )
        {
          unsigned row = order[ rowIndex ];
          unsigned column = getLeadingZeros( row );
          for ( unsigned subRowIndex = 0; subRowIndex < rowIndex; ++subRowIndex )
          {
            unsigned subRow = order[ subRowIndex ];
            rowOperation
            (
              subRow,
              row,
              -mat[ subRow ][ column ]
            );
          }
        }

      } // reducedRowEcholon

      //-------------------------------------------------------------
      // Return the determinant of the matrix.
      // Recursive function.
      //-------------------------------------------------------------
      TYPE determinant() const
      {
        TYPE result = static_cast< TYPE >( 0 );

        // Must have a square matrix to even bother.
        assert( rows == columns );

        if ( rows > 2 )
        {
          int sign = 1;
          for ( unsigned column = 0; column < columns; ++column )
          {
            TYPE subDeterminant;

            matrix subMatrix = mat( *this, 0, column );

            subDeterminant  = subMatrix.determinant();
            subDeterminant *= mat[ 0 ][ column ];

            if ( sign > 0 )
              result += subDeterminant;
            else
              result -= subDeterminant;

            sign = -sign;
          }
        }
        else
        {
          result = ( mat[ 0 ][ 0 ] * mat[ 1 ][ 1 ] )
                 - ( mat[ 0 ][ 1 ] * mat[ 1 ][ 0 ] );
        }

        return result;

      } // determinant

      //-------------------------------------------------------------
      // Calculate a dot product between this and an other matrix.
      //-------------------------------------------------------------
      TYPE dotProduct
      (
        matrix const & otherMatrix
      ) const
      {
        // Dimentions of each matrix must be the same.
        assert( rows == otherMatrix.rows );
        assert( columns == otherMatrix.columns );

        TYPE result = static_cast< TYPE >( 0 );
        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
          {
            result +=
              mat[ row ][ column ]
              * otherMatrix.mat[ row ][ column ];
          }

        return result;

      } // dotProduct

      //-------------------------------------------------------------
      // Return the transpose of the matrix.
      //-------------------------------------------------------------
      matrix const getTranspose() const
      {
        matrix result( columns, rows );

        // Transpose the matrix by filling the result's rows will
        // these columns, and vica versa.
        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            result.mat[ column ][ row ] = mat[ row ][ column ];

        return result;

      } // transpose

      //-------------------------------------------------------------
      // Transpose the matrix.
      //-------------------------------------------------------------
      void transpose()
      {
        *this = getTranspose();
      }

      //-------------------------------------------------------------
      // Return inverse matrix.
      //-------------------------------------------------------------
      matrix const getInverse() const
      {
        // Concatenate the identity matrix onto this matrix.
        matrix inverseMatrix
          (
            *this,
            IdentityMatrix< TYPE >( rows, columns ),
            TO_RIGHT
          );

        // Row reduce this matrix.  This will result in the identity
        // matrix on the left, and the inverse matrix on the right.
        inverseMatrix.reducedRowEcholon();

        // Copy the inverse matrix data back to this matrix.
        matrix result
        (
          inverseMatrix.getSubMatrix
          (
            0,
            rows - 1,
            columns,
            columns + columns - 1,
            inverseMatrix.order
          )
        );

        return result;

      } // invert


      //-------------------------------------------------------------
      // Invert this matrix.
      //-------------------------------------------------------------
      void invert()
      {
        *this = getInverse();

      } // invert

      //-----------------------------------------------------------------------------
      // *** ADDED: save function used to save a matrix class to a file (on disk).
      //-----------------------------------------------------------------------------
      void saveMatrix(std::string fname, matrix<TYPE> G, bool withSize){
      std::ofstream fout;
      fout.open(fname.c_str());
        if(withSize){
          fout << G.nr() << std::endl;
          }
        for ( unsigned int  i = 0 ; i < G.nr(); i++ ){
           for(unsigned int  j=0; j < G.nc(); j++){
      	fout << G(i,j) << " ";
           }
           fout << std::endl;
        }
      fout.close();
      }

      //-----------------------------------------------------------------------------
      // *** ADDED: max function over the matrix elements.
      //-----------------------------------------------------------------------------
      TYPE max(matrix<TYPE> mat){
      double xmax = -1e300;
        for ( unsigned int  i = 0 ; i < mat.nr(); i++ ){
           for(unsigned int  j=0; j < mat.nc(); j++){
      	if(xmax<mat(i,j))
      		xmax=mat(i,j);
           }
        }
      return xmax;
      }
      //-----------------------------------------------------------------------------
      // *** ADDED: min function over the matrix elements.
      //-----------------------------------------------------------------------------
      TYPE min(matrix<TYPE> mat){
      double xmin = 1e300;
        for ( unsigned int  i = 0 ; i < mat.nr(); i++ ){
           for(unsigned int  j=0; j < mat.nc(); j++){
      	if(xmin>mat(i,j))
      		xmin=mat(i,j);
           }
        }
      return xmin;
      }
      //---------------------------------------------------------------------------
      // *** ADDED: replace each element of the input matrix  by its logarithm (neperian).
      //---------------------------------------------------------------------------
      matrix<TYPE>  ComputeLogMatrix(matrix<TYPE>  GM){
      for (unsigned int  i = 0 ; i < GM.nr(); i++ ){
         for(unsigned int  j=i; j < GM.nc(); j++){
        	 GM(i,j)=log(GM(i,j)+1e-300);
        	 GM(j,i)=GM(i,j);
         }
      }
      return (GM);
      }


      //----------------------------------------------------------------------------------------
      // *** ADDED: replace each element of the input matrix  by its 'normalized' logarithm (neperian).
      //-----------------------------------------------------------------------------------------
      matrix<TYPE> ComputeLogNormalizedMatrix(matrix<TYPE> GM){
      for (unsigned int  i = 0 ; i < GM.nr(); i++ ){
         for(unsigned int  j=i; j < GM.nc(); j++){
      	   GM(i,j)=log(GM(i,j)+1e-300);
      	   GM(j,i)=GM(i,j);
         }
      }
      double logmax=GM.max();
      double logmin=GM.min();
      for (unsigned int  i = 0 ; i < GM.nr(); i++ ){
         for(unsigned int  j=i; j < GM.nc(); j++){
      	   GM(i,j)=exp((GM(i,j)-logmin)/(logmax-logmin));
      	   GM(j,i)=GM(i,j);
         }
      }
      return (GM);
      }

      //----------------------------------------------------------------------------------------------------
      // *** replace each element of the input matrix M by its 't' power, with t=1.0/(log(max(M)-log(min(M)).
      //----------------------------------------------------------------------------------------------------
      matrix<TYPE>  ComputePowerNormalizedMatrix( matrix<TYPE>  GM, double *t){
      double logmin = 1e300;
      double logmax = -1e300, x;
      for ( unsigned int  i = 0 ; i < GM.nr(); i++ ){
         for(unsigned int  j=i; j < GM.nc(); j++){
            x=log(GM(i,j)+1e-300);
            if (logmin>x)
            	  logmin=x;
            if (logmax<x)
            	  logmax=x;
         }
      }
      *t=1.0/(logmax-logmin);
      for ( unsigned int  i = 0 ; i < GM.nc(); i++ ){
         for(unsigned int  j=i; j < GM.nr(); j++){
        	 GM(i,j)=pow(GM(i,j),*t);
        	 GM(j,i)=GM(i,j);
         }
      }
      return (GM);
      }

      //=======================================================================
      // Operators.
      //=======================================================================

      //-------------------------------------------------------------
      // *** ADDED: acess (r,c) operator.
      //-------------------------------------------------------------
      // get value at (row,col)
       TYPE& operator()(int row, int col)
       {
           return mat[row][col];
       }

      //-------------------------------------------------------------
      // Add by an other matrix.
      //-------------------------------------------------------------
      matrix const operator +
      (
        matrix const & otherMatrix
      ) const
      {
        assert( otherMatrix.rows == rows );
        assert( otherMatrix.columns == columns );

        matrix result( rows, columns );

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            result.mat[ row ][ column ] =
              mat[ row ][ column ]
              + otherMatrix.mat[ row ][ column ];

        return result;
      }

      //-------------------------------------------------------------
      // Add self by an other matrix.
      //-------------------------------------------------------------
      matrix const & operator +=
      (
        matrix const & otherMatrix
      )
      {
        *this = *this + otherMatrix;
        return *this;
      }

      //-------------------------------------------------------------
      // Subtract by an other matrix.
      //-------------------------------------------------------------
      matrix const operator -
      (
        matrix const & otherMatrix
      ) const
      {
        assert( otherMatrix.rows == rows );
        assert( otherMatrix.columns == columns );

        matrix result( rows, columns );

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            result.mat[ row ][ column ] =
              mat[ row ][ column ]
              - otherMatrix.mat[ row ][ column ];

        return result;
      }

      //-------------------------------------------------------------
      // Subtract self by an other matrix.
      //-------------------------------------------------------------
      matrix const & operator -=
      (
        matrix const & otherMatrix
      )
      {
        *this = *this - otherMatrix;
        return *this;
      }

      //-------------------------------------------------------------
      // *** ADDED: Matrix element by element multiplication.
      //-------------------------------------------------------------
//      matrix const operator x
//      (
//        matrix const & otherMatrix
//      ) const
//      {
//        assert( otherMatrix.rows == rows );
//        assert( otherMatrix.columns == columns );
//
//        matrix result( rows, columns );
//
//        for ( unsigned row = 0; row < rows; ++row )
//          for ( unsigned column = 0; column < columns; ++column )
//          {
//            result.mat[ row ][ column ] = mat[row][column]*otherMatrix.mat[row][column];
//           }
//
//        return result;
//      }

      //-------------------------------------------------------------
      // Matrix multiplication.
      //-------------------------------------------------------------
      matrix const operator *
      (
        matrix const & otherMatrix
      ) const
      {
        TYPE const ZERO = static_cast< TYPE >( 0 );

        assert( otherMatrix.rows == columns );

        matrix result( rows, otherMatrix.columns );

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < otherMatrix.columns; ++column )
          {
            result.mat[ row ][ column ] = ZERO;

            for ( unsigned index = 0; index < columns; ++index )
              result.mat[ row ][ column ] +=
                mat[ row ][ index ]
                * otherMatrix.mat[ index ][ column ];
          }

        return result;
      }

      //-------------------------------------------------------------
      // Multiply self by matrix.
      //-------------------------------------------------------------
      matrix const & operator *=
      (
        matrix const & otherMatrix
      )
      {
        *this = *this * otherMatrix;
        return *this;
      }

      //-------------------------------------------------------------
      // Multiply by scalar constant.
      //-------------------------------------------------------------
      matrix const operator *
      (
        TYPE const & scalar
      ) const
      {
        matrix result( rows, columns );

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            result.mat[ row ][ column ] = mat[ row ][ column ] * scalar;

        return result;
      }

      //-------------------------------------------------------------
      // Multiply self by scalar constant.
      //-------------------------------------------------------------
      matrix const & operator *=
      (
        TYPE const & scalar
      )
      {
        *this = *this * scalar;
        return *this;
      }

      //-------------------------------------------------------------
      // Copy matrix.
      //-------------------------------------------------------------
      matrix & operator =
      (
        matrix const & otherMatrix
      )
      {
        if ( this == &otherMatrix )
          return *this;

        // Release memory currently in use.
        deallocate( rows, columns );

        rows    = otherMatrix.rows;
        columns = otherMatrix.columns;
        allocate( rows, columns );

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            mat[ row ][ column ] =
            otherMatrix.mat[ row ][ column ];

        return *this;
      }

      //-------------------------------------------------------------
      // Copy matrix data from array.
      // Although matrix data is two dimensional, this copy function
      // assumes the previous row is immediately followed by the next
      // row's data.
      //
      // Example for 3x2 matrix:
      //     int const data[ 3 * 2 ] =
      //     {
      //       1, 2, 3,
      //       4, 5, 6
      //     };
      //    Matrix< int > matrix( 3, 2 );
      //    matrix = data;
      //-------------------------------------------------------------
      matrix & operator =
      (
        TYPE const * data
      )
      {
        unsigned index = 0;

        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            mat[ row ][ column ] = data[ index++ ];

        return *this;
      }

      //-----------------------------------------------------------------------
      // Return true if this matrix is the same of parameter.
      //-----------------------------------------------------------------------
      bool operator ==
      (
        matrix const & value
      ) const
      {
        bool isEqual = true;
        for ( unsigned row = 0; row < rows; ++row )
          for ( unsigned column = 0; column < columns; ++column )
            if ( mat[ row ][ column ] != value.mat[ row ][ column ] )
              isEqual = false;

        return isEqual;
      }

      //-----------------------------------------------------------------------
      // Return true if this matrix is NOT the same of parameter.
      //-----------------------------------------------------------------------
      bool operator !=
      (
        matrix const & value
      ) const
      {
        return !( *this == value );
      }

      //-------------------------------------------------------------
      // Constructor for empty matrix.
      // Only useful if matrix is being assigned (i.e. copied) from
      // somewhere else sometime after construction.
      //-------------------------------------------------------------
      matrix()
      :
        rows( 0 ),
        columns( 0 )
      {
        allocate( 0, 0 );
      }

      //-------------------------------------------------------------
      // Constructor using rows and columns.
      //-------------------------------------------------------------
      matrix
      (
        unsigned rowsParameter,
        unsigned columnsParameter
      )
      :
        rows( rowsParameter ),
        columns( columnsParameter )
      {
        TYPE const ZERO = static_cast< TYPE >( 0 );

        // Allocate memory for new matrix.
        allocate( rows, columns );

        // Fill matrix with zero.
        for ( unsigned row = 0; row < rows; ++row )
        {
          order[ row ] = row;

          for ( unsigned column = 0; column < columns; ++column )
            mat[ row ][ column ] = ZERO;
        }
      }

      //-------------------------------------------------------------
      // This constructor will allow the creation of a matrix based off
      // an other matrix.  It can copy the matrix entirely, or omitted a
      // row/column.
      //-------------------------------------------------------------
      matrix
      (
        matrix const & copyMatrix,
        unsigned omittedRow    = INT_MAX,
        unsigned omittedColumn = INT_MAX
      )
      {
        // Start with the number of rows/columns from matrix to be copied.
        rows    = copyMatrix.getRows();
        columns = copyMatrix.getColumns();

        // If a row is omitted, then there is one less row.
        if ( INT_MAX != omittedRow  )
          rows--;

        // If a column is omitted, then there is one less column.
        if ( INT_MAX != omittedColumn )
          columns--;

        // Allocate memory for new matrix.
        allocate( rows, columns );

        unsigned rowIndex = 0;
        for ( unsigned row = 0; row < rows; ++row )
        {
          // If this row is to be skipped...
          if ( rowIndex == omittedRow )
            rowIndex++;

          // Set default order.
          order[ row ] = row;

          unsigned columnIndex = 0;
          for ( unsigned column = 0; column < columns; ++column )
          {
            // If this column is to be skipped...
            if ( columnIndex == omittedColumn )
              columnIndex++;

            mat[ row ][ column ] = copyMatrix.mat[ rowIndex ][ columnIndex ];

            columnIndex++;
          }

          ++rowIndex;
        }

      }

      //-------------------------------------------------------------
      // Constructor to concatenate two matrices.  Concatenation
      // can be done to the right, or to the bottom.
      //   A = [B | C]
      //-------------------------------------------------------------
      matrix
      (
        matrix const & copyMatrixA,
        matrix const & copyMatrixB,
        Position position = TO_RIGHT
      )
      {
        unsigned rowOffset    = 0;
        unsigned columnOffset = 0;

        if ( TO_RIGHT == position )
          columnOffset = copyMatrixA.columns;
        else
          rowOffset = copyMatrixA.rows;

        rows    = copyMatrixA.rows    + rowOffset;
        columns = copyMatrixA.columns + columnOffset;

        // Allocate memory for new matrix.
        allocate( rows, columns );

        for ( unsigned row = 0; row < copyMatrixA.rows; ++row )
          for ( unsigned column = 0; column < copyMatrixA.columns; ++column )
            mat[ row ][ column ] = copyMatrixA.mat[ row ][ column ];

        for ( unsigned row = 0; row < copyMatrixB.rows; ++row )
          for ( unsigned column = 0; column < copyMatrixB.columns; ++column )
            mat[ row + rowOffset ][ column + columnOffset ] =
              copyMatrixB.mat[ row ][ column ];
      }

      //-------------------------------------------------------------
      // Destructor.
      //-------------------------------------------------------------
      ~matrix()
      {
        // Release memory.
        deallocate( rows, columns );
      }

  };

//=============================================================================
// Class for identity matrix.
//=============================================================================
template< class TYPE >
  class IdentityMatrix : public matrix< TYPE >
  {
    public:
      IdentityMatrix
      (
        unsigned rowsParameter,
        unsigned columnsParameter
      )
      :
        matrix< TYPE >( rowsParameter, columnsParameter )
      {
        TYPE const ZERO = static_cast< TYPE >( 0 );
        TYPE const ONE  = static_cast< TYPE >( 1 );

        for ( unsigned row = 0; row < matrix< TYPE >::rows; ++row )
        {
          for ( unsigned column = 0; column < matrix< TYPE >::columns; ++column )
            if ( row == column )
            	matrix< TYPE >::mat[ row ][ column ] = ONE;
            else
            	matrix< TYPE >::mat[ row ][ column ] = ZERO;
        }
      }
  };

//-----------------------------------------------------------------------------
// Stream operator used to convert matrix class to a string.
//-----------------------------------------------------------------------------
template< class TYPE >
  std::ostream & operator<<
  (
    // Stream data to place string.
    std::ostream & stream,

    // A matrix.
	matrix< TYPE > const & M
  )
  {
    for ( unsigned row = 0; row < M.getRows(); ++row )
    {
      for ( unsigned column = 0; column < M.getColumns(); ++column )
        stream << "\t" << M.get( row , column );

      stream << std::endl;
    }

    return stream;
  }



#endif // _MATRIX_H_
