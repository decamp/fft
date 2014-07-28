/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

/**
 * Performs a Fast Discrete Cosine Transform on an array of real values.
 * <p>
 * Not thread safe.
 *
 * @author Philip DeCamp
 */
public class FastCosineTransform {

    private final int mDim;
    private final int mBits;

    private final double[] mWeight;
    private final double[] mInvWeight;
    private final double[] mWorkA;


    /**
     * The sole argument, <tt>dim</tt> indicates the size
     * of vectors on which the transform will operate.  This
     * must be a power-of-two.
     * <p>
     * Memory footprint is a bit over 48 * dim bytes.
     *
     * @param dim Size of vector on which this transform operates.  Must be power-of-two.
     * @throws IllegalArgumentException if dim is smaller than 2, too large, or not a power-of-two.
     */
    public FastCosineTransform( int dim ) {
        mDim = dim;
        mBits = FastFourierTransform.computeBitNum( dim );

        mWeight = new double[dim * 2];
        mInvWeight = new double[dim * 2];
        mWorkA = new double[dim * 2];

        computeWeightVectors( dim, mWeight, mInvWeight );
    }

    /**
     * Performs a Fast Discrete Cosine Transform on an array of real values.
     * <p>
     * Not thread safe.
     *
     * @param a       Input array of real values with size [dim].
     * @param aOff    Offset into array <tt>a</tt>
     * @param inverse Set to <tt>true</tt> to perform inverse transform.
     * @param out     Output matrix where DCT coeffs are stored.  Must have space for <tt>dim</tt> values.
     * @param outOff  Offset into array <tt>out</tt>
     */
    public void apply( double[] a, int aOff, boolean inverse, double[] out, int outOff ) {
        if( !inverse ) {
            shuffle1( a, aOff, mDim, mBits, mWorkA );
            FastFourierTransform.transform( mWorkA, 0, mDim, false );
            shuffle2( mWorkA, mWeight, mDim, out, outOff );
        } else {
            invShuffle1( a, aOff, mInvWeight, mDim, mBits, mWorkA );
            FastFourierTransform.transform( mWorkA, 0, mDim, true );
            invShuffle2( mWorkA, mDim, out, outOff );
        }
    }



    private static void computeWeightVectors( int dim, double[] out, double[] outInv ) {
        final double s = 1.0 / dim;
        out[0] = 1.0;
        out[1] = 0.0;
        outInv[0] = s;
        outInv[1] = 0.0;

        for( int i = 1; i < dim; i++ ) {
            double cos = Math.cos(  i * Math.PI * 0.5 / dim );
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
            int rowIn = ( FastFourierTransform.reverse( rowOut2 ) >>> shift );

            //Butterfly shuffle.
            rowIn = rowIn + (rowIn / dim) * (dim2 - 2 * rowIn - 1);
            out[rowOut2] = a[rowIn + offA];
            out[rowOut2 + 1] = 0.0;
        }
    }

    /**
     * 1. Apply weight vector       <br>
     * 2. Drop imaginary components <br>
     */
    private static void shuffle2( double[] a, double[] w, int dim, double[] out, int offOut ) {
        final int dim2 = 2 * dim;

        for( int y2 = 0; y2 < dim2; y2 += 2 ) {
            out[offOut++] = a[y2] * w[y2] - a[y2 + 1] * w[y2 + 1];
        }
    }

    /**
     * 1. Apply weights       <br>
     * 2. Reverse-bit shuffle <br>
     */
    private static void invShuffle1( double[] a, int offA, double[] w, int dim, int bits, double[] out ) {
        final int shift = 31 - bits;

        for( int rowIn = 0; rowIn < dim; rowIn++ ) {
            int rowOut2 = ( FastFourierTransform.reverse( rowIn ) >>> shift );
            double v = a[offA + rowIn];

            out[rowOut2] = v * w[rowIn * 2];
            out[rowOut2 + 1] = v * w[rowIn * 2 + 1];
        }
    }

    /**
     * 1. Drop imaginary component   <br>
     * 2. Reverse butterfly shuffle  <br>
     */
    private static void invShuffle2( double[] a, int dim, double[] out, int offOut ) {
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Reverse butterfly shuffle.
            final int rowO = i2 + (i2 / dim) * (dim2 - 1 - 2 * i2);
            out[rowO + offOut] = a[i2];
        }
    }

}