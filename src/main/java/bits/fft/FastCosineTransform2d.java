/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

/**
 * Performs a Fast Discrete Cosine Transform on an square matrix of real values.
 * <p>
 * Not thread safe.
 *
 * @author Philip DeCamp
 */
public class FastCosineTransform2d {

    private final int mDim;
    private final int mBits;

    private final double[] mWeight;
    private final double[] mInvWeight;
    private final double[] mWorkA;
    private final double[] mWorkB;


    /**
     * The sole argument, <tt>dim</tt>, indicates the size
     * of the matrices on which this transform will operate,
     * which must be a power-of-two. For example, if dim == 8,
     * then all input arrays must be have a length of at least 8*8.
     * <p>
     * Memory footprint is a bit over 32 * ( dim * dim + dim ) bytes.
     *
     * @param dim Size of one-side of square matrix on which this transform operates.  Must be power-of-two.
     * @throws IllegalArgumentException if dim is smaller than 2, too large, or not a power-of-two.
     */
    public FastCosineTransform2d( int dim ) {
        mDim  = dim;
        mBits = FastFourierTransform2d.computeBitNum( dim );

        mWeight    = new double[dim * 2];
        mInvWeight = new double[dim * 2];
        mWorkA     = new double[dim * dim * 2];
        mWorkB     = new double[dim * dim * 2];

        computeWeightVectors( dim, mWeight, mInvWeight );
    }

    /**
     * Performs a Fast Discrete Cosine Transform on a square matrix of real values.
     * <p>
     * Samples must be stored in a tight-packed format. For a 2x2 matrix: <br>
     * <tt>[... r_0_0, r_1_0, r_0_1, r_1_1, ...] </tt><br>,
     * where <tt>r_m_n</tt> is the real-valued element at position [m,n].
     *
     * @param a       Input matrix of real values with size [dim,dim].
     * @param aOff    Offset into array <tt>a</tt>
     * @param inverse Set to <tt>true</tt> to perform inverse transform.
     * @param out     Output matrix where DCT coeffs are stored.  Must have space for <tt>dim*dim</tt> values.
     * @param outOff  Offset into array<tt>out</tt>
     */
    public void apply( double[] a, int aOff, boolean inverse, double[] out, int outOff ) {
        if( !inverse ) {
            shuffle1( a, aOff, mDim, mBits, mWorkA );
            FastFourierTransform2d.transform( mWorkA, 0, mDim, false );
            shuffle2( mWorkA, mWeight, mDim, mBits, mWorkB );
            FastFourierTransform2d.transform( mWorkB, 0, mDim, false );
            shuffle3( mWorkB, mWeight, mDim, out, outOff );
        } else {
            invShuffle1( a, aOff, mInvWeight, mDim, mBits, mWorkA );
            FastFourierTransform2d.transform( mWorkA, 0, mDim, true );
            invShuffle2( mWorkA, mInvWeight, mDim, mBits, mWorkB );
            FastFourierTransform2d.transform( mWorkB, 0, mDim, true );
            invShuffle3( mWorkB, mDim, out, outOff );
        }
    }



    private static void computeWeightVectors( int dim, double[] out, double[] outInv ) {
        final double s = 1.0 / dim;
        out[0] = 1.0;
        out[1] = 0.0;
        outInv[0] = s;
        outInv[1] = 0.0;

        for( int i = 1; i < dim; i++ ) {
            double cos = Math.cos(  i * Math.PI * 0.5 * s );
            double sin = Math.sqrt( 1.0 - cos * cos );

            out[i*2  ] =  2.0 * cos;
            out[i*2+1] = -2.0 * sin;

            outInv[i*2  ] = s * cos;
            outInv[i*2+1] = s * sin;
        }
    }

    /**
     * 1. Drop complex components <br>
     * 2. Butterfly shuffle       <br>
     * 3. Reverse-bit shuffle     <br>
     */
    private static void shuffle1( double[] a, int offA, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int rowOut2 = 0; rowOut2 < dim2; rowOut2 += 2 ) {
            //Bit reversal shuffle.
            int rowIn = (FastFourierTransform2d.reverse( rowOut2 ) >>> shift);

            //Butterfly shuffle.  
            //Equivalent to: if(x < dim / 2) {ia = x * 2;} else {ia = dim*2 - x*2 - 1;}
            //Example where dim = 8: [0 1 2 3 4 5 6 7] -> [0 2 4 6 7 5 3 1] 
            rowIn = rowIn + (rowIn / dim) * (dim2 - 2 * rowIn - 1);

            //Translate row A position into array element.
            int indA = rowIn + offA;

            for( int y = 0; y < dim; y++ ) {
                int vec = y * dim2;
                out[rowOut2 + vec] = a[indA + y * dim];
                out[rowOut2 + vec + 1] = 0.0;
            }
        }
    }

    /**
     * 1. Apply weight vector       <br>
     * 2. Drop imaginary components <br>
     * 3. Transpose                 <br>
     * 4. Butterfly shuffle         <br>
     * 5. Reverse-bit shuffle       <br>
     */
    private static void shuffle2( double[] a, double[] w, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int rowOut2 = 0; rowOut2 < dim2; rowOut2 += 2 ) {
            //Bit reversal shuffle.
            int rowA = (FastFourierTransform2d.reverse( rowOut2 ) >>> shift);

            //Butterfly shuffle.  
            //Equivalent to: if(x < dim / 2) {ia = x * 2;} else {ia = dim*2 - x*2 - 1;}
            //Example where dim = 8: [0 1 2 3 4 5 6 7] -> [0 2 4 6 7 5 3 1] 
            rowA = rowA + (rowA / dim) * (dim2 - 2 * rowA - 1);

            //Translate row A position into array offsets.
            //Note that A is being transposed.
            int indA = rowA * dim2;

            //Iterater through each column.
            for( int col2 = 0; col2 < dim2; col2 += 2 ) {
                //Perform complex multiplaction; throw away imaginary component.
                out[rowOut2 + col2 * dim] = a[indA + col2] * w[col2] - a[indA + col2 + 1] * w[col2 + 1];
                out[rowOut2 + col2 * dim + 1] = 0.0;
            }
        }
    }

    /**
     * 1. Apply weight vector.       <br>
     * 2. Transpose                  <br>
     * 3. Drop imaginary components. <br>
     *
     * @param a      Input: complex matrix of size [dim,dim]
     * @param w      Weight vector of size [dim] needed to compute DCT from FFT.
     * @param dim    Dimension of transform.  Must be power-of-two.
     * @param out    Output: real-valued matrix of size [dim,dim].
     * @param offOut Offset into <tt>out</tt>.
     */
    private static void shuffle3( double[] a, double[] w, int dim, double[] out, int offOut ) {
        final int dim2 = 2 * dim;

        for( int y2 = 0; y2 < dim2; y2 += 2 ) {
            double vReal = w[y2];
            double vImag = w[y2 + 1];

            for( int x = 0; x < dim; x++ ) {
                int ii = y2 + x * dim2;
                out[offOut++] = a[ii] * vReal - a[ii + 1] * vImag;
            }
        }
    }

    /**
     * 1. transform weights,  <br>
     * 2. reverse-bit shuffle <br>
     */
    private static void invShuffle1( double[] a, int offA, double[] w, int dim, int bits, double[] out ) {
        final int dim2 = dim * 2;
        final int shift = 31 - bits;

        for( int rowIn = 0; rowIn < dim; rowIn++ ) {
            int indA = offA + rowIn;
            int indOut = (FastFourierTransform2d.reverse( rowIn ) >>> shift);
            double wReal = w[rowIn * 2];
            double wImag = w[rowIn * 2 + 1];

            for( int col = 0; col < dim; col++ ) {
                double v = a[indA + col * dim];
                out[indOut + col * dim2] = wReal * v;
                out[indOut + col * dim2 + 1] = wImag * v;
            }
        }
    }

    /**
     * 1. Drop imaginary component  <br>
     * 2. Reverse butterfly shuffle <br>
     * 3. Transpose                 <br>
     * 4. Transform weights         <br>
     * 5. Reverse-bit shuffle       <br>
     */
    private static void invShuffle2( double[] a, double[] w, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Bit reverse shuffle.
            final int col2 = FastFourierTransform2d.reverse( i2 ) >>> shift;
            final int indA = col2 * dim;

            final double wReal = w[col2];
            final double wImag = w[col2 + 1];

            for( int j2 = 0; j2 < dim2; j2 += 2 ) {
                //Reverse butterfly shuffle.
                int rowO = j2 + (j2 / dim) * (dim2 - 1 - 2 * j2);
                double val = a[j2 + indA];

                out[i2 + rowO * dim2] = val * wReal;
                out[i2 + rowO * dim2 + 1] = val * wImag;
            }
        }
    }

    /**
     * 1. Drop imaginary component   <br>
     * 2. Reverse butterfly shuffle  <br>
     */
    private static void invShuffle3( double[] a, int dim, double[] out, int offOut ) {
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Reverse butterfly shuffle.
            final int rowO = i2 + (i2 / dim) * (dim2 - 1 - 2 * i2);
            final int indO = offOut + rowO * dim;

            for( int j = 0; j < dim; j++ ) {
                out[j + indO] = a[i2 + j * dim2];
            }
        }
    }

}